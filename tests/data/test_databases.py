import os

from plinder.data import databases


def test_make_db_limits_external_tool_threads(tmp_path, monkeypatch):
    commands: list[list[str]] = []
    monkeypatch.setattr(databases, "run", commands.append)

    databases.make_db(
        input_dir=tmp_path / "input",
        output_dir=tmp_path / "output",
        db="foldseek",
        threads=3,
    )

    assert len(commands) == 2
    for command in commands:
        index = command.index("--threads")
        assert command[index + 1] == "3"
    foldseek_createdb = commands[0]
    assert foldseek_createdb[foldseek_createdb.index("--coord-store-mode") + 1] == "2"


def _touch_database(database, *, foldseek=False):
    database.parent.mkdir(exist_ok=True, parents=True)
    database.with_suffix(".dbtype").touch()
    database.with_suffix(".index").touch()
    if foldseek:
        database.parent.joinpath(f"{database.name}_ss.dbtype").touch()
        database.parent.joinpath(f"{database.name}_ca.dbtype").touch()


def test_exact_foldseek_target_uses_native_cluster_search_db(tmp_path, monkeypatch):
    full = tmp_path / "subdbs" / "holo_foldseek" / "holo_foldseek"
    _touch_database(full, foldseek=True)
    commands: list[list[str]] = []

    def fake_run(command, **_kwargs):
        commands.append(command)
        if command[1] == "linclust":
            _touch_database(full.parent / "exact_clusters")
        elif command[1] == "createclusearchdb":
            _touch_database(full.parent / "clustered", foldseek=True)
        elif command[1] == "createindex":
            (full.parent / "clustered.idx.dbtype").touch()

    monkeypatch.setattr(databases, "run", fake_run)

    report = databases.make_exact_search_db(
        full_db=full,
        aln_type="foldseek",
        tmp_dir=tmp_path / "scratch",
        threads=4,
    )

    assert [command[1] for command in commands] == [
        "linclust",
        "createclusearchdb",
        "createindex",
    ]
    linclust = commands[0]
    assert linclust[linclust.index("--min-seq-id") + 1] == "1.0"
    assert linclust[linclust.index("-c") + 1] == "1.0"
    cluster_db = commands[1]
    assert cluster_db[cluster_db.index("--compressed") + 1] == "0"
    assert report["search_target"] == "clustered"
    assert report["conversion_target"] == "clustered"
    assert report["source_index"]["name"] == "holo_foldseek.index"
    assert report["portable"] is True
    assert report["compressed_search_target"] is False


def test_exact_mmseqs_target_stores_expansion_alignments(tmp_path, monkeypatch):
    full = tmp_path / "subdbs" / "holo_mmseqs" / "holo_mmseqs"
    _touch_database(full)
    commands: list[list[str]] = []

    def fake_run(command, **_kwargs):
        commands.append(command)
        output_by_command = {
            "linclust": "exact_clusters",
            "createsubdb": "representatives",
            "align": "cluster_alignments",
        }
        if command[1] in output_by_command:
            _touch_database(full.parent / output_by_command[command[1]])
        elif command[1] == "createindex":
            (full.parent / "representatives.idx.dbtype").touch()

    monkeypatch.setattr(databases, "run", fake_run)

    report = databases.make_exact_search_db(
        full_db=full,
        aln_type="mmseqs",
        tmp_dir=tmp_path / "scratch",
        threads=2,
    )

    assert [command[1] for command in commands] == [
        "linclust",
        "createsubdb",
        "align",
        "createindex",
    ]
    cluster_align = commands[2]
    assert "-a" in cluster_align
    assert "--min-seq-id" not in cluster_align
    assert "-c" not in cluster_align
    assert cluster_align[cluster_align.index("--add-self-matches") + 1] == "1"
    createsubdb = commands[1]
    assert createsubdb[createsubdb.index("--subdb-mode") + 1] == "0"
    assert report["cluster_alignments"] == "cluster_alignments"


def test_portable_foldseek_subdb_copies_each_database(tmp_path, monkeypatch):
    full = tmp_path / "full" / "foldseek"
    full.parent.mkdir()
    full.with_suffix(".lookup").write_text("0 1abc_A 0\n")
    output = tmp_path / "selected"
    output.mkdir()
    commands: list[list[str]] = []
    monkeypatch.setattr(databases, "run", commands.append)

    missing = databases.make_sub_db(
        {"1abc_A"},
        full,
        output,
        "foldseek",
        portable=True,
    )

    assert missing == []
    assert len(commands) == 1
    assert all(
        command[command.index("--subdb-mode") + 1] == "0" for command in commands
    )


def test_make_sub_db_selects_only_first_model_per_chain(tmp_path, monkeypatch):
    full = tmp_path / "full" / "foldseek"
    full.parent.mkdir()
    full.with_suffix(".lookup").write_text(
        "0 pdb_00001abc_xyz-enrich_A 0\n"
        "1 pdb_00001abc_xyz-enrich_A 1\n"
        "2 pdb_00001abc_xyz-enrich_B 1\n"
    )
    output = tmp_path / "selected"
    output.mkdir()
    monkeypatch.setattr(databases, "run", lambda *_args, **_kwargs: None)

    missing = databases.make_sub_db(
        {
            "pdb_00001abc_xyz-enrich_A",
            "pdb_00001abc_xyz-enrich_B",
        },
        full,
        output,
        "foldseek",
        portable=True,
    )

    assert missing == []
    assert (output / "ids.tsv").read_text() == (
        "0 pdb_00001abc_xyz-enrich_A 0\n2 pdb_00001abc_xyz-enrich_B 1\n"
    )


def test_database_link_normalization_makes_directory_self_contained(tmp_path):
    root = tmp_path / "database"
    root.mkdir()
    internal_target = root / "member"
    internal_target.write_text("member")
    external_target = tmp_path / "parent-header"
    external_target.write_text("header")
    internal_link = root / "internal-link"
    internal_link.symlink_to(internal_target)
    external_link = root / "external-link"
    external_link.symlink_to(external_target)

    report = databases._make_database_directory_portable(root)

    assert report == {"external_links_copied": 1, "internal_links_relativized": 1}
    assert internal_link.is_symlink()
    assert not internal_link.readlink().is_absolute()
    assert external_link.is_file() and not external_link.is_symlink()
    assert external_link.read_text() == "header"
    assert not databases._has_external_database_links(root)


def test_database_link_checks_skip_runtime_alignment_trees(tmp_path):
    root = tmp_path / "database"
    runtime = root / "aln"
    runtime.mkdir(parents=True)
    external_target = tmp_path / "completed-alignment.parquet"
    external_target.write_text("alignment")
    runtime_link = runtime / "1abc.parquet"
    runtime_link.symlink_to(external_target)

    report = databases._make_database_directory_portable(root)

    assert report == {"external_links_copied": 0, "internal_links_relativized": 0}
    assert runtime_link.is_symlink()
    assert not databases._has_external_database_links(root)


def test_database_identifiers_follow_selected_index_keys(tmp_path):
    database = tmp_path / "selected"
    database.with_suffix(".lookup").write_text(
        "0\t1abc_A\t0\n1\t2def_B\t0\n2\t3ghi_C\t0\n"
    )
    database.with_suffix(".index").write_text("0\t0\t10\n2\t10\t20\n")

    assert databases.database_identifiers(database) == {"1abc_A", "3ghi_C"}


def test_install_database_directory_replaces_partial_shared_output(tmp_path):
    source = tmp_path / "scratch" / "foldseek"
    source.mkdir(parents=True)
    (source / "foldseek.dbtype").write_text("complete")
    source_payload = tmp_path / "scratch" / "payload"
    source_payload.write_text("portable")
    (source / "external-link").symlink_to(source_payload)
    target = tmp_path / "shared" / "foldseek"
    target.mkdir(parents=True)
    (target / "partial").write_text("stale")

    databases.install_database_directory(source, target)

    assert (target / "foldseek.dbtype").read_text() == "complete"
    assert (target / "external-link").read_text() == "portable"
    assert not (target / "external-link").is_symlink()
    assert not (target / "partial").exists()
    assert not (target.parent / ".foldseek.installing").exists()
    assert not (target.parent / ".foldseek.previous").exists()


def test_install_database_directory_preserves_generated_results(tmp_path):
    source = tmp_path / "scratch" / "foldseek"
    source.mkdir(parents=True)
    (source / "foldseek.dbtype").write_text("new database")
    target = tmp_path / "shared" / "foldseek"
    alignments = target / "aln"
    alignments.mkdir(parents=True)
    (alignments / "1abc.parquet").write_text("completed search")
    (target / "foldseek.dbtype").write_text("old database")

    databases.install_database_directory(
        source,
        target,
        preserve_directories=("aln", "mapped_aln"),
    )

    assert (target / "foldseek.dbtype").read_text() == "new database"
    assert (target / "aln" / "1abc.parquet").read_text() == "completed search"
    assert not (target.parent / ".foldseek.previous").exists()


def test_make_sub_dbs_builds_complete_backend_in_scratch(tmp_path, monkeypatch):
    full_db = tmp_path / "full" / "foldseek"
    full_db.parent.mkdir()
    full_db.with_suffix(".lookup").write_text("0 pdb_00001abc_xyz-enrich_A 0\n")
    db_dir = tmp_path / "shared" / "subdbs"
    db_dir.mkdir(parents=True)
    original_missing = tmp_path / "original-missing.json"
    original_missing.write_text('{"old": true}\n')
    os.link(original_missing, db_dir / "missing.json")
    scratch = tmp_path / "scratch"
    build_paths = []

    def fake_make_sub_db(_ids, _full_db, subdb, _aln_type, *, portable):
        build_paths.append(subdb)
        assert portable is True
        database = subdb / subdb.name
        _touch_database(database, foldseek=True)
        return []

    def fake_make_exact_search_db(*, full_db, **_kwargs):
        (full_db.parent / "exact_cluster.json").write_text("{}")
        return {"status": "complete"}

    monkeypatch.setattr(databases, "make_sub_db", fake_make_sub_db)
    monkeypatch.setattr(databases, "make_exact_search_db", fake_make_exact_search_db)
    monkeypatch.setattr(
        databases, "_completed_exact_search_manifest", lambda *_args: None
    )

    databases.make_sub_dbs(
        db_dir,
        {"holo_foldseek": full_db},
        identifiers_by_database={"holo_foldseek": {"pdb_00001abc_xyz-enrich_A"}},
        tmp_dir=scratch,
        threads=4,
    )

    assert build_paths == [scratch / "holo_foldseek/build/holo_foldseek"]
    target = db_dir / "holo_foldseek"
    assert (target / "holo_foldseek.dbtype").is_file()
    assert (target / "selection.json").is_file()
    assert (target / "exact_cluster.json").is_file()
    assert original_missing.read_text() == '{"old": true}\n'
    assert not os.path.samefile(original_missing, db_dir / "missing.json")
