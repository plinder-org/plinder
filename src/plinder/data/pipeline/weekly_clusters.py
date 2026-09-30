"""Extend fixed representative covers with changed release entries."""

from __future__ import annotations

from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Any

import duckdb
import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

from plinder.core.release import RELEASE_PATHS
from plinder.core.scores.metrics import CHEMICAL_CLUSTER_SUMMARY_COLUMNS
from plinder.data.databases import run
from plinder.data.protein_clusters import SEQUENCE_CLUSTER_SCHEMA


def _score_rows(
    paths: list[Path], query_ids: pd.Series, columns: list[str], *, numeric: bool
) -> pd.DataFrame:
    if not paths or query_ids.empty:
        return pd.DataFrame(columns=columns)
    ids = pd.DataFrame({"id": query_ids.drop_duplicates()})
    sources = ", ".join(
        f"'{str(path).replace(chr(39), chr(39) * 2)}'" for path in paths
    )
    selected = ", ".join(f"scores.{column}" for column in columns)
    query_column = columns[0]
    with duckdb.connect() as connection:
        connection.register("changed_queries", ids)
        rows = connection.execute(
            f"SELECT {selected} FROM read_parquet([{sources}]) AS scores "
            f"JOIN changed_queries ON scores.{query_column} = changed_queries.id"
        ).df()
    if numeric:
        rows[query_column] = rows[query_column].astype("Int64")
        rows[columns[1]] = rows[columns[1]].astype("Int64")
    return rows


def _best_representatives(
    scores: pd.DataFrame,
    representatives: pd.DataFrame,
    *,
    threshold: int,
    directed: bool,
    previous_labels: pd.Series | None = None,
) -> pd.DataFrame:
    """Prefer a primary-threshold hit, then the directed-cover 50% fallback."""
    if scores.empty or representatives.empty:
        return pd.DataFrame(columns=["query_node", "label", "centroid"])
    candidates = scores.merge(
        representatives,
        left_on="target_node",
        right_on="centroid",
        how="inner",
    )
    floor = min(threshold, 50) if directed else threshold
    candidates = candidates.loc[candidates["similarity"].ge(floor)].copy()
    if candidates.empty:
        return pd.DataFrame(columns=["query_node", "label", "centroid"])
    candidates["primary"] = candidates["similarity"].ge(threshold)
    candidates["previous"] = (
        candidates["label"].eq(candidates["query_node"].map(previous_labels))
        if previous_labels is not None
        else False
    )
    return (
        candidates.sort_values(
            ["query_node", "previous", "primary", "similarity", "centroid"],
            ascending=[True, False, False, False, True],
        )
        .drop_duplicates("query_node")[["query_node", "label", "centroid"]]
        .reset_index(drop=True)
    )


def extend_ligand_clusters(
    *, base: Path, data_dir: Path, affected: set[str]
) -> pd.DataFrame:
    """Keep existing covers fixed and assign changed ligands to their centroids."""
    annotation = pd.read_parquet(
        data_dir / "index/annotation_table.parquet",
        columns=[
            "entry_pdb_id",
            "ligand_id",
            "ligand_smiles",
            "ligand_is_proper",
            "system_type",
        ],
    )
    chemistry = pd.read_parquet(
        data_dir / "fingerprints/ligand_similarity_annotations.parquet",
        columns=["ligand_rdkit_canonical_smiles", "ligand_smiles_id"],
    ).rename(columns={"ligand_rdkit_canonical_smiles": "ligand_smiles"})
    annotation = annotation.merge(
        chemistry, on="ligand_smiles", how="left", validate="many_to_one"
    )
    current = annotation.drop_duplicates("ligand_id").set_index("ligand_id")
    changed = current.loc[current["entry_pdb_id"].isin(affected)]
    proper = changed["ligand_is_proper"].fillna(False) & changed["system_type"].eq(
        "holo"
    )
    changed_ids = changed.index[proper].astype(str)
    previous = pd.read_parquet(base / "index/ligand_clusters.parquet").set_index(
        "ligand_id"
    )
    result = previous.reindex(current.index.astype(str)).copy()
    result.loc[changed.index.astype(str), :] = pd.NA

    if len(changed_ids):
        shards = {ligand_id[1:3] for ligand_id in changed_ids}
        score_paths = [
            data_dir / "exports/ligand_similarity_scores" / f"{shard}.parquet"
            for shard in sorted(shards)
        ]
        score_paths = [path for path in score_paths if path.is_file()]
        direct = _score_rows(
            score_paths,
            pd.Series(changed_ids, dtype="string"),
            [
                "query_ligand_id",
                "target_ligand_id",
                "pocket_qcov",
                "pli_qcov",
                "sucos_shape",
            ],
            numeric=False,
        )
        direct["sucos_shape_pocket_qcov"] = np.floor(
            direct["sucos_shape"] * direct["pocket_qcov"] / 100
        )
        smiles = current["ligand_smiles_id"].dropna().astype("int64")
        changed_smiles = smiles.loc[smiles.index.intersection(changed_ids)]
        changed_lookup = changed_smiles.rename_axis("ligand_id").reset_index(
            name="smiles_id"
        )
        chemical: dict[str, pd.DataFrame] = {}
        for metric, directory in (
            ("tanimoto_similarity_ecfp4_1024", "ligand_scores"),
            ("jaccard_similarity_mhfp6_2048", "mhfp6_scores"),
        ):
            chemical[metric] = _score_rows(
                sorted((data_dir / directory).glob("*.parquet")),
                changed_smiles.reset_index(drop=True),
                ["query_ligand_id", "target_ligand_id", metric],
                numeric=True,
            )

        for column in previous.columns:
            if not column.endswith(("__set_cover", "__directed_set_cover")):
                continue
            metric, threshold_text, _, mode = column.split("__")
            threshold = int(threshold_text)
            centroid_column = f"{column}__is_centroid"
            centroids = previous.loc[
                previous[centroid_column].fillna(False), [column]
            ].reset_index()
            centroids.columns = ["centroid", "label"]
            centroids = centroids.loc[centroids["centroid"].isin(current.index)]
            if mode == "directed_set_cover":
                scores = direct[["query_ligand_id", "target_ligand_id", metric]].rename(
                    columns={
                        "query_ligand_id": "query_node",
                        "target_ligand_id": "target_node",
                        metric: "similarity",
                    }
                )
            else:
                centroid_smiles = centroids.merge(
                    smiles.rename("smiles_id"),
                    left_on="centroid",
                    right_index=True,
                    how="inner",
                )
                edges = chemical[metric].merge(
                    changed_lookup,
                    left_on="query_ligand_id",
                    right_on="smiles_id",
                    how="inner",
                )
                scores = edges.merge(
                    centroid_smiles,
                    left_on="target_ligand_id",
                    right_on="smiles_id",
                    how="inner",
                    suffixes=("", "_centroid"),
                )[["ligand_id", "centroid", metric]].rename(
                    columns={
                        "ligand_id": "query_node",
                        "centroid": "target_node",
                        metric: "similarity",
                    }
                )
            matches = _best_representatives(
                scores,
                centroids,
                threshold=threshold,
                directed=mode == "directed_set_cover",
                previous_labels=previous[column],
            ).set_index("query_node")
            centroid_labels = centroids.set_index("centroid")["label"]
            for ligand_id in changed_ids:
                # A revised cover anchor keeps its ID until the scheduled
                # full rebuild; the reference cover is fixed for this period.
                if ligand_id in centroid_labels.index:
                    label = centroid_labels.at[ligand_id]
                    representative = ligand_id
                elif ligand_id in matches.index:
                    label = matches.at[ligand_id, "label"]
                    representative = matches.at[ligand_id, "centroid"]
                else:
                    label = f"weekly:{ligand_id}"
                    representative = ligand_id
                result.at[ligand_id, column] = label
                result.at[ligand_id, centroid_column] = representative == ligand_id
                if mode == "directed_set_cover":
                    result.at[ligand_id, f"{column}__coverage_count"] = pd.NA
                    result.at[ligand_id, f"{column}__coverage_fraction"] = pd.NA

    for metric, summary in CHEMICAL_CLUSTER_SUMMARY_COLUMNS.items():
        if summary not in result:
            continue
        cover = f"{metric}__90__ligand__set_cover"
        result.loc[changed_ids, summary] = result.loc[changed_ids, cover]
        label_pdb = annotation.loc[
            annotation["ligand_is_proper"].fillna(False)
            & annotation["system_type"].eq("holo"),
            ["entry_pdb_id", "ligand_id"],
        ].merge(result[[summary]].reset_index(), on="ligand_id", how="left")
        counts = label_pdb.groupby(summary)["entry_pdb_id"].nunique()
        result[f"{summary}_num_pdb_ids"] = result[summary].map(counts).astype("Int32")
    return result.reset_index().rename(columns={"index": "ligand_id"})


def extend_interface_clusters(
    *, base: Path, data_dir: Path, affected: set[str]
) -> pd.DataFrame:
    """Assign changed interfaces to the fixed whole- and half-interface covers."""
    annotation = pd.read_parquet(
        data_dir / "index/interface_annotation_table.parquet",
        columns=["entry_pdb_id", "system_id"],
    ).drop_duplicates("system_id")
    current_ids = pd.Index(annotation["system_id"].astype(str), name="system_id")
    changed_ids = annotation.loc[
        annotation["entry_pdb_id"].isin(affected), "system_id"
    ].astype(str)
    previous = pd.read_parquet(base / "index/interface_clusters.parquet").set_index(
        "system_id"
    )
    result = previous.reindex(current_ids).copy()
    result.loc[changed_ids, :] = pd.NA
    if changed_ids.empty:
        return result.reset_index()

    membership = pd.read_parquet(
        data_dir / "index/interface_membership.parquet",
        columns=[
            "system_id",
            "representative_system_id",
            "side_1_half_interface_id",
            "side_2_half_interface_id",
        ],
    ).set_index("system_id")
    changed_membership = membership.loc[changed_ids]
    changed_nodes = pd.Series(
        pd.unique(changed_membership.to_numpy().ravel()), dtype="string"
    ).dropna()
    paths = [
        data_dir / "interface_scores" / f"shard={shard}.parquet"
        for shard in sorted({node[1:3] for node in changed_nodes})
    ]
    paths = [path for path in paths if path.is_file()]
    edges = _score_rows(
        paths,
        changed_nodes,
        ["query_system", "target_system", "metric", "similarity"],
        numeric=False,
    )
    current_nodes = set(membership.to_numpy().ravel())
    for column in previous.columns:
        if not column.endswith("_directed_set_cover"):
            continue
        metric, threshold_text, kind = column.split("__")
        threshold = int(threshold_text)
        source = (
            base
            / "interface_sampling/directed_set_cover"
            / f"metric={metric}"
            / f"threshold={threshold}.parquet"
        )
        labels = (
            pq.ParquetFile(source)
            .read(columns=["system_id", "centroid_system_id", "label"])
            .to_pandas()
        )
        centroids = labels.loc[
            labels["system_id"].eq(labels["centroid_system_id"]),
            ["system_id", "label"],
        ].rename(columns={"system_id": "centroid"})
        centroids = centroids.loc[centroids["centroid"].isin(current_nodes)]
        centroid_labels = centroids.set_index("centroid")["label"]
        score_rows = edges.loc[edges["metric"].eq(metric)].rename(
            columns={
                "query_system": "query_node",
                "target_system": "target_node",
            }
        )[["query_node", "target_node", "similarity"]]
        membership_column = (
            "representative_system_id"
            if kind == "directed_set_cover"
            else "side_1_half_interface_id"
            if kind == "chain_1_directed_set_cover"
            else "side_2_half_interface_id"
        )
        previous_labels = pd.Series(
            previous[column].reindex(changed_membership.index).to_numpy(),
            index=changed_membership[membership_column].to_numpy(),
        ).dropna()
        previous_labels = previous_labels.loc[
            ~previous_labels.index.duplicated(keep="first")
        ]
        matches = _best_representatives(
            score_rows,
            centroids,
            threshold=threshold,
            directed=True,
            previous_labels=previous_labels,
        ).set_index("query_node")
        for system_id, member in changed_membership.iterrows():
            node = member[membership_column]
            if node in centroid_labels.index:
                label = centroid_labels.at[node]
            elif node in matches.index:
                label = matches.at[node, "label"]
            else:
                label = f"weekly:{node}"
            result.at[system_id, column] = label
    return result.reset_index()


def extend_protein_clusters(
    *,
    base: Path,
    data_dir: Path,
    affected: set[str],
    backend: str,
    threshold: float,
    coverage: float,
    threads: int,
    scratch_dir: Path,
) -> Path:
    """Search changed chains against fixed protein-cluster representatives."""
    if backend not in {"mmseqs", "foldseek"}:
        raise ValueError(f"unknown protein clustering backend: {backend}")
    artifact = (
        "protein_sequence_clusters"
        if backend == "mmseqs"
        else "protein_structure_clusters"
    )
    old = pd.read_parquet(base / RELEASE_PATHS[artifact])
    chains = pd.read_parquet(
        data_dir / "index/entry_chains.parquet",
        columns=["entry_pdb_id", "chain_asym_id", "chain_auth_id", "chain_sequence"],
        filters=[("chain_receptor_type", "==", "protein")],
    )
    key = ["entry_pdb_id", "chain_asym_id"]
    current = chains[key].merge(old, on=key, how="left", validate="one_to_one")
    current["is_representative"] = current["is_representative"].astype("boolean")
    changed = current["entry_pdb_id"].isin(affected)
    current.loc[
        changed,
        [
            "representative_entry_pdb_id",
            "representative_chain_asym_id",
            "is_representative",
        ],
    ] = pd.NA
    current.loc[changed, "status"] = (
        "missing_sequence" if backend == "mmseqs" else "not_in_foldseek_db"
    )
    query_chains = chains.loc[
        chains["entry_pdb_id"].isin(affected) & chains["chain_auth_id"].notna()
    ].copy()
    representatives = old.loc[
        old["is_representative"].astype("boolean").fillna(False), key
    ]
    representatives = representatives.merge(chains, on=key, how="inner")
    representatives = representatives.loc[representatives["chain_auth_id"].notna()]
    if backend == "mmseqs":
        query_chains = query_chains.loc[query_chains["chain_sequence"].notna()]
        representatives = representatives.loc[representatives["chain_sequence"].notna()]

    def chain_name(row: Any) -> str:
        if backend == "mmseqs":
            return f"{row.entry_pdb_id}_{row.chain_auth_id}"
        return f"pdb_0000{row.entry_pdb_id}_xyz-enrich_{row.chain_auth_id}"

    query_names = {
        chain_name(row): (row.entry_pdb_id, row.chain_asym_id)
        for row in query_chains.itertuples(index=False)
    }
    target_names = {
        chain_name(row): (row.entry_pdb_id, row.chain_asym_id)
        for row in representatives.itertuples(index=False)
    }
    target_by_key = {chain_key: name for name, chain_key in target_names.items()}
    old_by_key = old.set_index(key)[
        ["representative_entry_pdb_id", "representative_chain_asym_id"]
    ]
    preferred = {
        name: target_by_key.get(tuple(old_by_key.loc[chain_key]))
        for name, chain_key in query_names.items()
        if chain_key in old_by_key.index
    }
    base_db = data_dir / "dbs" / backend / backend
    overlay_db = data_dir / "dbs/weekly_delta" / backend / backend
    output = data_dir / RELEASE_PATHS[artifact]
    lookup: dict[str, str] = {}
    matched: dict[str, str] = {}
    if query_names and target_names:
        scratch_dir.mkdir(parents=True, exist_ok=True)
        with TemporaryDirectory(prefix=f"weekly-{backend}-", dir=scratch_dir) as temp:
            work = Path(temp)
            query_db = (
                overlay_db if overlay_db.with_suffix(".dbtype").is_file() else base_db
            )
            with Path(f"{query_db}.lookup").open() as handle:
                for line in handle:
                    db_id, name, _ = line.rstrip("\n").split("\t", 2)
                    if name in query_names and name not in lookup:
                        lookup[name] = db_id
            if lookup:
                keys = work / "query.keys"
                keys.write_text("".join(f"{db_id}\n" for db_id in lookup.values()))
                run(
                    [
                        backend,
                        "createsubdb",
                        str(keys),
                        str(query_db),
                        str(work / "query"),
                        "--subdb-mode",
                        "0",
                    ]
                )
                shadowed_path = data_dir / "manifests/weekly_shadowed_entries.parquet"
                shadowed = (
                    set(pd.read_parquet(shadowed_path, columns=["pdb_id"]).pdb_id)
                    if shadowed_path.is_file()
                    else set()
                )
                fields = "query,target,fident,qcov,tcov"
                if backend == "foldseek":
                    fields += ",lddt"
                parts = []
                for label, target_db, names in (
                    (
                        "base",
                        base_db,
                        {
                            name
                            for name, chain in target_names.items()
                            if chain[0] not in shadowed
                        },
                    ),
                    (
                        "overlay",
                        overlay_db,
                        {
                            name
                            for name, chain in target_names.items()
                            if chain[0] in shadowed
                        },
                    ),
                ):
                    if not names or not target_db.with_suffix(".dbtype").is_file():
                        continue
                    target_lookup: dict[str, str] = {}
                    with Path(f"{target_db}.lookup").open() as handle:
                        for line in handle:
                            db_id, name, _ = line.rstrip("\n").split("\t", 2)
                            if name in names and name not in target_lookup:
                                target_lookup[name] = db_id
                    if not target_lookup:
                        continue
                    target_keys = work / f"{label}.keys"
                    target_keys.write_text(
                        "".join(f"{db_id}\n" for db_id in target_lookup.values())
                    )
                    target = work / f"target-{label}"
                    run(
                        [
                            backend,
                            "createsubdb",
                            str(target_keys),
                            str(target_db),
                            str(target),
                            "--subdb-mode",
                            "0",
                        ]
                    )
                    hits_db = work / f"hits-{label}"
                    command = [
                        backend,
                        "search",
                        str(work / "query"),
                        str(target),
                        str(hits_db),
                        str(work / f"tmp-{label}"),
                        "-c",
                        str(coverage),
                        "--cov-mode",
                        "0",
                        "--max-seqs",
                        "10000",
                        "--threads",
                        str(threads),
                    ]
                    command += (
                        ["--min-seq-id", str(threshold)]
                        if backend == "mmseqs"
                        else [
                            "--lddt-threshold",
                            str(threshold),
                            "--alignment-type",
                            "2",
                            "-a",
                            "1",
                        ]
                    )
                    run(command)
                    hits_tsv = work / f"hits-{label}.tsv"
                    run(
                        [
                            backend,
                            "convertalis",
                            str(work / "query"),
                            str(target),
                            str(hits_db),
                            str(hits_tsv),
                            "--format-output",
                            fields,
                            "--threads",
                            str(threads),
                        ]
                    )
                    if hits_tsv.stat().st_size:
                        parts.append(
                            pd.read_csv(hits_tsv, sep="\t", names=fields.split(","))
                        )
                hits = (
                    pd.concat(parts, ignore_index=True)
                    if parts
                    else pd.DataFrame(columns=fields.split(","))
                )
                quality = "fident" if backend == "mmseqs" else "lddt"
                hits = hits.loc[
                    hits[quality].ge(threshold)
                    & hits["qcov"].ge(coverage)
                    & hits["tcov"].ge(coverage)
                    & hits["target"].isin(target_names)
                ]
                # Keep the prior assignment when it still passes. This makes
                # revisions stable even if another representative scores better.
                hits["preferred"] = hits["target"].eq(hits["query"].map(preferred))
                best = hits.sort_values(
                    ["query", "preferred", quality, "target"],
                    ascending=[True, False, False, True],
                ).drop_duplicates("query")
                matched = dict(zip(best["query"], best["target"]))
    indexed = current.set_index(key)
    for name, chain_key in query_names.items():
        if name not in lookup:
            continue
        representative = target_names.get(matched.get(name, ""), chain_key)
        indexed.loc[chain_key, "representative_entry_pdb_id"] = representative[0]
        indexed.loc[chain_key, "representative_chain_asym_id"] = representative[1]
        indexed.loc[chain_key, "is_representative"] = representative == chain_key
        indexed.loc[chain_key, "status"] = "weekly_assigned"
    table = pa.Table.from_pandas(
        indexed.reset_index()[SEQUENCE_CLUSTER_SCHEMA.names],
        schema=SEQUENCE_CLUSTER_SCHEMA,
        preserve_index=False,
    )
    output.parent.mkdir(parents=True, exist_ok=True)
    staged = output.with_suffix(".weekly.parquet")
    pq.write_table(table, staged, compression="zstd")
    staged.replace(output)
    return output
