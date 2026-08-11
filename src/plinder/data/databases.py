# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
import hashlib
import json
import os
import shutil
import subprocess as sp
from collections.abc import Iterator
from errno import EACCES, EPERM, EXDEV
from pathlib import Path
from typing import TYPE_CHECKING, Any, Dict, Optional

import pandas as pd

from plinder.core.utils.log import setup_logger

if TYPE_CHECKING:
    from plinder.core.scores.entries import EntryView

LOG = setup_logger(__name__)

EXACT_CLUSTER_IDENTITY = 1.0
EXACT_CLUSTER_COVERAGE = 1.0
DATABASE_RUNTIME_DIRECTORIES = {"aln", "mapped_aln"}


def run(cmd: list[str], *, cwd: Path | None = None) -> None:
    LOG.info(" ".join(cmd))
    sp.check_call(cmd, stderr=sp.STDOUT, cwd=cwd)


def make_db(
    input_dir: Path,
    output_dir: Path,
    db: str,
    tmp_dir: Path = Path("tmp"),
    threads: int = 1,
) -> None:
    """
    Create full databases. input_dir is a misnomer
    for mmseqs because it wants the full path to the
    pdb_seqres.txt.gz or pdb_uniprot.fasta. output_dir is a misnomer in both
    cases because it wants the full path to the filename
    of the database, not just the directory it writes to.

    Parameters
    ----------
    input_dir : Path
        location of input files for database
        (foldseek: input_dir = adir [/ **/*-enrich.cif.gz] for apo/holo and [/AF-*-F1-model_v4.cif] for pred)
        (mmseqs: input_dir = adir / seqres / pdb_seqres.txt.gz for apo/holo and / uniprot / pdb_uniprot.txt.gz for pred)
    output_dir : Path
        location of full final database (including file name)
    db : str
        name of db in ["foldseek", "mmseqs"]
    tmp_dir : Path, default=Path("tmp")
        scratch directory
    threads : int, default=1
        Maximum number of threads used by Foldseek or MMseqs.
    """
    create_db(input_dir, output_dir, db, threads=threads)
    create_db_index(output_dir, db, tmp_dir=tmp_dir, threads=threads)


def create_db(
    input_dir: Path,
    output_dir: Path,
    db: str,
    *,
    threads: int = 1,
) -> None:
    """Create a Foldseek/MMseqs database without indexing it."""
    full_db = str(output_dir / db)
    output_dir.mkdir(exist_ok=True, parents=True)
    _remove_database_prefix(Path(full_db))
    cmd = [
        db,
        "createdb",
        str(input_dir),
        full_db,
        "--threads",
        str(threads),
    ]
    if db == "foldseek":
        # Foldseek only needs one coordinate per residue.  Make the compact
        # C-alpha/delta representation explicit so a future default change
        # cannot inflate the build or portable search database.
        cmd.extend(["--chain-name-mode", "1", "--coord-store-mode", "2"])
    run(cmd)


def create_db_index(
    output_dir: Path,
    db: str,
    *,
    tmp_dir: Path = Path("tmp"),
    threads: int = 1,
) -> None:
    """Index an existing Foldseek/MMseqs database."""
    full_db = str(output_dir / db)
    run(
        [
            db,
            "createindex",
            full_db,
            str(tmp_dir) + f"_{db}",
            "--threads",
            str(threads),
        ]
    )


def _database_files_exist(database: Path, aln_type: str) -> bool:
    """Return whether the minimum files for a sequence/structure DB exist."""
    required = [database.with_suffix(".dbtype")]
    if aln_type == "foldseek":
        required.extend(
            [
                Path(f"{database}_ss.dbtype"),
                Path(f"{database}_ca.dbtype"),
            ]
        )
    return all(path.is_file() for path in required)


def created_database_is_complete(database: Path, aln_type: str) -> bool:
    """Return whether ``createdb`` finished all lookup and data outputs."""
    return _database_files_exist(database, aln_type) and all(
        path.is_file()
        for path in [database.with_suffix(".index"), database.with_suffix(".lookup")]
    )


def install_database_directory(
    source: Path,
    target: Path,
    *,
    preserve_directories: tuple[str, ...] = (),
) -> None:
    """Copy a completed scratch DB and swap it into shared storage.

    Generated search results can share a backend directory with its target
    database.  ``preserve_directories`` moves those directories through the
    atomic swap instead of copying them to scratch or deleting them with the
    retired database generation.
    """
    staging = target.parent / f".{target.name}.installing"
    backup = target.parent / f".{target.name}.previous"
    target.parent.mkdir(exist_ok=True, parents=True)
    if backup.exists() and not target.exists():
        backup.rename(target)
    for path in [staging, backup]:
        if path.exists():
            shutil.rmtree(path)
    shutil.copytree(source, staging, symlinks=True)
    _make_database_directory_portable(staging)
    if target.exists():
        target.rename(backup)
    try:
        for name in preserve_directories:
            preserved = backup / name
            if not preserved.exists():
                continue
            destination = staging / name
            if destination.exists():
                raise FileExistsError(
                    f"cannot preserve {preserved}: {destination} already exists"
                )
            preserved.rename(destination)
        staging.rename(target)
    except BaseException:
        for name in preserve_directories:
            preserved = staging / name
            destination = backup / name
            if preserved.exists() and backup.exists() and not destination.exists():
                preserved.rename(destination)
        if backup.exists() and not target.exists():
            backup.rename(target)
        raise
    if backup.exists():
        shutil.rmtree(backup)


def _remove_database_prefix(database: Path) -> None:
    """Remove one generated MMseqs/Foldseek DB and all of its side files."""
    for path in database.parent.glob(f"{database.name}*"):
        if path.is_dir() and not path.is_symlink():
            shutil.rmtree(path)
        else:
            path.unlink(missing_ok=True)


def _database_bundle_paths(root: Path) -> Iterator[Path]:
    """Iterate database assets without walking large runtime result trees."""
    for current, directories, filenames in os.walk(root, topdown=True):
        current_path = Path(current)
        if current_path == root:
            directories[:] = [
                name
                for name in directories
                if name not in DATABASE_RUNTIME_DIRECTORIES
            ]
        yield from (current_path / name for name in directories)
        yield from (current_path / name for name in filenames)


def _make_database_directory_portable(root: Path) -> dict[str, int]:
    """Replace external DB links with copies and internal links with relatives."""
    root_resolved = root.resolve()
    copied = 0
    relativized = 0
    for path in sorted(_database_bundle_paths(root)):
        if not path.is_symlink():
            continue
        target = path.resolve(strict=True)
        path.unlink()
        if target.is_relative_to(root_resolved):
            path.symlink_to(
                os.path.relpath(target, start=path.parent),
                target_is_directory=target.is_dir(),
            )
            relativized += 1
        elif target.is_dir():
            shutil.copytree(target, path)
            copied += 1
        else:
            shutil.copy2(target, path)
            copied += 1
    return {"external_links_copied": copied, "internal_links_relativized": relativized}


def _has_external_database_links(root: Path) -> bool:
    root_resolved = root.resolve()
    for path in _database_bundle_paths(root):
        if not path.is_symlink():
            continue
        if path.readlink().is_absolute():
            return True
        try:
            target = path.resolve(strict=True)
        except FileNotFoundError:
            return True
        if not target.is_relative_to(root_resolved):
            return True
    return False


def _link_or_copy_database_file(source: Path, destination: Path) -> None:
    """Install one database file without duplicating it on the ingest volume."""
    if source.is_symlink():
        destination.symlink_to(source.readlink())
        return
    try:
        os.link(source, destination)
    except OSError as exc:
        if exc.errno not in {EACCES, EPERM, EXDEV}:
            raise
        shutil.copy2(source, destination)


def publish_search_database_bundle(
    *,
    source_root: Path,
    target_root: Path,
    aln_type: str,
) -> dict[str, int | str]:
    """Atomically publish the portable files needed for custom searches.

    The full Foldseek source database and exact-clustering build artifacts are
    not needed after ``createclusearchdb``. MMseqs retains its full target and
    representative-to-member alignments because ``expandaln`` followed by
    member-level ``align`` requires both.
    """
    if aln_type not in {"foldseek", "mmseqs"}:
        raise ValueError(f"unsupported alignment type: {aln_type}")
    full_db = source_root / source_root.name
    manifest = _completed_exact_search_manifest(full_db, aln_type)
    if manifest is None:
        raise ValueError(f"incomplete exact-cluster database: {source_root}")

    prefixes = {
        str(manifest["search_target"]),
        str(manifest["conversion_target"]),
    }
    if aln_type == "mmseqs":
        prefixes.add(str(manifest["cluster_alignments"]))
    sources = {source_root / "exact_cluster.json"}
    for prefix in prefixes:
        sources.update(
            path
            for path in source_root.glob(f"{prefix}*")
            if path.is_file() or path.is_symlink()
        )
    if _has_external_database_links(source_root):
        raise ValueError(f"non-portable database links in {source_root}")

    staging = target_root.parent / f".{target_root.name}.installing"
    backup = target_root.parent / f".{target_root.name}.previous"
    target_root.parent.mkdir(exist_ok=True, parents=True)
    for path in (staging, backup):
        if path.exists():
            shutil.rmtree(path)
    staging.mkdir()
    try:
        for source in sorted(sources):
            _link_or_copy_database_file(source, staging / source.name)
        if _has_external_database_links(staging):
            raise ValueError(f"published database bundle is not portable: {staging}")
        if target_root.exists():
            target_root.rename(backup)
        staging.rename(target_root)
    except BaseException:
        if staging.exists():
            shutil.rmtree(staging)
        if backup.exists() and not target_root.exists():
            backup.rename(target_root)
        raise
    if backup.exists():
        shutil.rmtree(backup)

    regular_files = [path for path in target_root.iterdir() if path.is_file()]
    return {
        "alignment_type": aln_type,
        "file_count": len(regular_files),
        "apparent_size": sum(path.stat().st_size for path in regular_files),
        "search_target": str(manifest["search_target"]),
        "conversion_target": str(manifest["conversion_target"]),
    }


def _path_signature(path: Path, *, portable: bool = False) -> dict[str, int | str]:
    stat = path.stat()
    signature: dict[str, int | str] = {
        "size": stat.st_size,
        "mtime_ns": stat.st_mtime_ns,
    }
    signature["name" if portable else "path"] = (
        path.name if portable else str(path.resolve())
    )
    return signature


def _identifier_digest(identifiers: set[str]) -> str:
    digest = hashlib.sha256()
    for identifier in sorted(identifiers):
        digest.update(identifier.encode())
        digest.update(b"\n")
    return digest.hexdigest()


def _database_entry_count(database: Path) -> int:
    with database.with_suffix(".index").open("rb") as handle:
        return sum(1 for _ in handle)


def database_identifiers(database: Path) -> set[str]:
    """Return identifiers for records actually selected into a database."""
    lookup = {}
    with database.with_suffix(".lookup").open() as handle:
        for line in handle:
            fields = line.split()
            if len(fields) >= 2:
                lookup[fields[0]] = fields[1]
    selected = set()
    with database.with_suffix(".index").open() as handle:
        for line in handle:
            fields = line.split()
            if fields and fields[0] in lookup:
                selected.add(lookup[fields[0]])
    return selected


def _read_json(path: Path) -> dict[str, Any] | None:
    try:
        payload = json.loads(path.read_text())
    except (OSError, ValueError, TypeError):
        return None
    return payload if isinstance(payload, dict) else None


def _write_json(path: Path, payload: dict[str, Any]) -> None:
    temporary = path.with_suffix(".tmp.json")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    temporary.replace(path)


def exact_search_database_paths(
    db_dir: Path,
    search_db: str,
    aln_type: str,
) -> tuple[Path, Path, Path | None]:
    """Return full target, representative target and expansion alignment DB."""
    root = db_dir / f"{search_db}_{aln_type}"
    full_target = root / root.name
    if aln_type == "foldseek":
        return full_target, root / "clustered", None
    if aln_type == "mmseqs":
        return full_target, root / "representatives", root / "cluster_alignments"
    raise ValueError(f"unsupported alignment type: {aln_type}")


def _completed_exact_search_manifest(
    full_db: Path, aln_type: str
) -> dict[str, Any] | None:
    """Return a valid exact-cluster manifest without mutating its database."""
    root = full_db.parent
    manifest = _read_json(root / "exact_cluster.json")
    try:
        source_signature = _path_signature(full_db.with_suffix(".index"), portable=True)
    except OSError:
        return None
    expected = {
        "alignment_type": aln_type,
        "identity": EXACT_CLUSTER_IDENTITY,
        "coverage": EXACT_CLUSTER_COVERAGE,
        "coverage_mode": 0,
        "source_index": source_signature,
        "portable": True,
        "compressed_search_target": False,
    }
    search_target = root / (
        "clustered" if aln_type == "foldseek" else "representatives"
    )
    expected_outputs = [
        search_target.with_suffix(".dbtype"),
        Path(f"{search_target}.idx.dbtype"),
    ]
    if aln_type == "mmseqs":
        expected_outputs.append((root / "cluster_alignments").with_suffix(".dbtype"))
    if (
        manifest is not None
        and all(manifest.get(key) == value for key, value in expected.items())
        and all(path.is_file() for path in expected_outputs)
        and not _has_external_database_links(root)
    ):
        return manifest
    return None


def make_exact_search_db(
    *,
    full_db: Path,
    aln_type: str,
    tmp_dir: Path,
    threads: int,
) -> dict[str, Any]:
    """Build and index an exactly deduplicated, expandable search target.

    The search cap is applied to representatives. Foldseek's native cluster
    search realigns every expanded structure. MMseqs uses a representative
    search followed by ``expandaln`` and a final member-level ``align``.
    """
    if aln_type not in {"foldseek", "mmseqs"}:
        raise ValueError(f"unsupported alignment type: {aln_type}")
    if not _database_files_exist(full_db, aln_type):
        raise FileNotFoundError(f"missing full {aln_type} target DB: {full_db}")

    root = full_db.parent
    manifest_path = root / "exact_cluster.json"
    source_signature = _path_signature(full_db.with_suffix(".index"), portable=True)
    expected = {
        "alignment_type": aln_type,
        "identity": EXACT_CLUSTER_IDENTITY,
        "coverage": EXACT_CLUSTER_COVERAGE,
        "coverage_mode": 0,
        "source_index": source_signature,
        "portable": True,
        "compressed_search_target": False,
    }
    search_target = root / (
        "clustered" if aln_type == "foldseek" else "representatives"
    )
    expansion_db = root / "cluster_alignments" if aln_type == "mmseqs" else None
    manifest = _completed_exact_search_manifest(full_db, aln_type)
    if manifest is not None:
        return manifest

    clusters = root / "exact_clusters"
    representatives = root / "representatives"
    cluster_alignments = root / "cluster_alignments"
    clustered = root / "clustered"
    for database in [clusters, representatives, cluster_alignments, clustered]:
        _remove_database_prefix(database)
    manifest_path.unlink(missing_ok=True)
    tmp_dir.mkdir(exist_ok=True, parents=True)

    run(
        [
            aln_type,
            "linclust",
            str(full_db),
            str(clusters),
            str(tmp_dir / "linclust"),
            "--min-seq-id",
            str(EXACT_CLUSTER_IDENTITY),
            "-c",
            str(EXACT_CLUSTER_COVERAGE),
            "--cov-mode",
            "0",
            "--threads",
            str(threads),
        ]
    )
    if aln_type == "foldseek":
        # Foldseek 10.941's cluster-search expansion segfaults in
        # mergeresultsbyset for compressed createclusearchdb targets.
        run(
            [
                "foldseek",
                "createclusearchdb",
                str(full_db),
                str(clusters),
                str(clustered),
                "--threads",
                str(threads),
                "--compressed",
                "0",
            ]
        )
        search_target = clustered
    else:
        run(
            [
                "mmseqs",
                "createsubdb",
                str(clusters),
                str(full_db),
                str(representatives),
                "--subdb-mode",
                "0",
            ]
        )
        # expandaln needs representative-to-member alignments, not only the
        # key-only cluster membership emitted by linclust.  Do not repeat the
        # exact-clustering thresholds here: linclust can cluster sequences
        # with terminal X residues that MMseqs cannot align at 100% coverage.
        # This database is only an expansion map; the post-expansion align
        # below applies the configured search e-value and coverage filters.
        run(
            [
                "mmseqs",
                "align",
                str(representatives),
                str(full_db),
                str(clusters),
                str(cluster_alignments),
                "-a",
                "-e",
                "1e100",
                "--add-self-matches",
                "1",
                "--threads",
                str(threads),
            ]
        )
        search_target = representatives

    run(
        [
            aln_type,
            "createindex",
            search_target.name,
            str((tmp_dir / "index").resolve()),
            "--threads",
            str(threads),
        ],
        cwd=root,
    )
    link_report = _make_database_directory_portable(root)
    member_count = _database_entry_count(full_db)
    representative_count = _database_entry_count(search_target)
    payload = {
        **expected,
        "search_target": search_target.name,
        "conversion_target": (
            search_target.name if aln_type == "foldseek" else full_db.name
        ),
        "cluster_database": clusters.name,
        "member_count": member_count,
        "representative_count": representative_count,
        "redundant_member_count": member_count - representative_count,
        "link_normalization": link_report,
    }
    if expansion_db is not None:
        payload["cluster_alignments"] = expansion_db.name
    _write_json(manifest_path, payload)
    return payload


def make_sub_db(
    entry_chain_ids: set[str],
    full_db: Path,
    sub_db: Path,
    aln_type: str,
    *,
    portable: bool = False,
) -> list[str]:
    """
    Requires a list of already generated annotation entries.
    Create reduced indexes from the full databases based
    on the passed entry_chain_ids.

    Parameters
    ----------
    entry_chain_ids : set[str]
        output from get_db_ids
    full_db : Path
        path to full database
    sub_db : Path
        path to sub database
    aln_type : str
        database name
    """
    # Map chain ids to full database ids
    sub_db_lookup_file = sub_db / "ids.tsv"
    database = sub_db / sub_db.name
    _remove_database_prefix(database)
    found = set()
    with open(full_db.parent / f"{full_db.name}.lookup", "r") as f:
        with open(sub_db_lookup_file, "w") as fw:
            for line in f:
                fields = line.split()
                # Foldseek emits one record per model but uses the same chain
                # identifier for every model.  PLINDER annotations use the
                # first model, and the lookup is ordered by model, so retain
                # only the first record for each requested chain.
                if fields[1] in entry_chain_ids and fields[1] not in found:
                    fw.write(line)
                    found.add(fields[1])
    # Create subdb
    if aln_type == "foldseek":
        run(
            [
                "foldseek",
                "createsubdb",
                str(sub_db_lookup_file),
                str(full_db),
                str(database),
                "--subdb-mode",
                "0" if portable else "1",
            ]
        )
    elif aln_type == "mmseqs":
        run(
            [
                "mmseqs",
                "createsubdb",
                str(sub_db_lookup_file),
                str(full_db),
                str(database),
            ]
        )

    missing = set(entry_chain_ids) - found
    return list(missing)


def get_ids_in_db(data_dir: Path, search_db: str, aln_type: str) -> pd.DataFrame:
    sub_db_file = data_dir / "dbs" / "subdbs" / f"{search_db}_{aln_type}" / "ids.tsv"
    return pd.read_csv(sub_db_file, sep="\t", header=None)


def get_db_ids(
    entries: Dict[str, "EntryView"],
    search_db: str,
    aln_type: str,
    entry_ids: Optional[list[str]] = None,
) -> set[str]:
    """
    get chain IDs for a given database from
    the provided entries.

    Parameters
    ----------
    entries : Dict[str, EntryView]
        the collection of entries
    search_db : str
        search_db in ["apo", "holo", "pred"]
    aln_type : str
        aln_type in ["foldseek", "mmseqs"]

    Returns
    -------
    db_ids : set[str]
        db IDs in the given database
    """
    db_ids = set()
    if entry_ids is None:
        entry_ids = list(entries.keys())
    for entry_id in entry_ids:
        for chain in entries[entry_id].chains_for_alignment(search_db, aln_type):
            db_ids.add(chain)
    return db_ids


def make_sub_dbs(
    db_dir: Path,
    db_sources: Dict[str, Path],
    entries: Dict[str, "EntryView"] | None = None,
    *,
    identifiers_by_database: dict[str, set[str]] | None = None,
    tmp_dir: Path | None = None,
    threads: int = 1,
) -> None:
    """
    Create the apo/holo subdbs for score
    generation.

    Parameters
    ----------
    db_dir : Path
        directory where subdbs get written
    db_sources : Dict[str, Path]
        map of database name to path to full database
    entries : Dict[str, EntryView], optional
        Map of entries used to derive identifiers for legacy apo/pred workflows.
    identifiers_by_database : dict[str, set[str]], optional
        Precomputed identifiers keyed by ``{search_db}_{alignment_type}``.
        V3 holo scoring uses this to select protein receptor chains directly
        from the normalized chain index without materializing all systems.
    """

    db_dir.mkdir(exist_ok=True, parents=True)
    report = {}
    cluster_report = {}
    for search_db_aln_type, full_db in db_sources.items():
        search_db, aln_type = search_db_aln_type.split("_")
        subdb = db_dir / search_db_aln_type
        subdb.mkdir(exist_ok=True)
        if (
            identifiers_by_database is not None
            and search_db_aln_type in identifiers_by_database
        ):
            requested_ids = identifiers_by_database[search_db_aln_type]
        elif entries is not None:
            requested_ids = get_db_ids(entries, search_db, aln_type)
        else:
            raise ValueError(f"no identifiers supplied for {search_db_aln_type}")
        selection_manifest = subdb / "selection.json"
        source_lookup = full_db.parent / f"{full_db.name}.lookup"
        selection = {
            "alignment_type": aln_type,
            "portable": True,
            "identifier_count": len(requested_ids),
            "identifier_sha256": _identifier_digest(requested_ids),
            "source_lookup": _path_signature(source_lookup),
        }
        cached = _read_json(selection_manifest)
        selected_database = subdb / subdb.name
        cluster_manifest = _completed_exact_search_manifest(selected_database, aln_type)
        complete = (
            cached is not None
            and all(cached.get(key) == value for key, value in selection.items())
            and _database_files_exist(selected_database, aln_type)
            and cluster_manifest is not None
        )
        if complete:
            assert cached is not None
            missing = list(cached.get("missing", []))
        else:
            backend_tmp = (tmp_dir or db_dir / "scratch") / search_db_aln_type
            working_subdb = subdb
            if tmp_dir is not None:
                working_subdb = backend_tmp / "build" / search_db_aln_type
                if working_subdb.exists():
                    shutil.rmtree(working_subdb)
                working_subdb.mkdir(exist_ok=True, parents=True)
            missing = sorted(
                make_sub_db(
                    requested_ids,
                    full_db,
                    working_subdb,
                    aln_type,
                    portable=True,
                )
            )
            working_selection = working_subdb / "selection.json"
            _write_json(working_selection, {**selection, "missing": missing})
            working_database = working_subdb / working_subdb.name
            cluster_manifest = make_exact_search_db(
                full_db=working_database,
                aln_type=aln_type,
                tmp_dir=backend_tmp / "work",
                threads=threads,
            )
            if tmp_dir is not None:
                install_database_directory(
                    working_subdb,
                    subdb,
                    preserve_directories=("aln", "mapped_aln"),
                )
        report[search_db_aln_type] = missing
        if cluster_manifest is None:
            raise RuntimeError(
                f"exact clustering did not complete for {search_db_aln_type}"
            )
        cluster_report[search_db_aln_type] = cluster_manifest
    with (db_dir / "missing.json").open("w") as f:
        json.dump(report, f)
    _write_json(db_dir / "exact_clusters.json", cluster_report)
