# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import TypeAlias

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.dataset as ds
import pyarrow.parquet as pq

from plinder.core.utils.schemas import STRUCTURE_LINK_SCHEMA

TableInput: TypeAlias = pd.DataFrame | str | Path

REQUIRED_SCORE_METRICS = (
    "pocket_fident",
    "protein_fident_weighted_sum",
    "protein_fident_qcov_weighted_sum",
    "protein_lddt_weighted_sum",
)
SCORE_COLUMNS = (
    "query_system",
    "query_ligand_id",
    "target_system",
    "metric",
    "similarity",
)
ANNOTATION_COLUMNS = ("system_id", "ligand_id", "ligand_is_proper")
HOLO_CHAIN_ANNOTATION_COLUMNS = (
    "entry_pdb_id",
    "ligand_is_proper",
    "ligand_protein_chains_asym_id",
)
ENTRY_CHAIN_COLUMNS = (
    "entry_pdb_id",
    "chain_asym_id",
    "chain_auth_id",
    "chain_entity_id",
    "chain_receptor_type",
    "chain_is_ligand_like",
)
BIOUNIT_CHAIN_COLUMNS = (
    "entry_pdb_id",
    "biounit_id",
    "chain_instance",
    "chain_asym_id",
    "chain_role",
    "chain_num_contacting_ions",
    "chain_num_contacting_artifacts",
    "chain_num_contacting_other_ligands",
)
ENTRY_METADATA_COLUMNS = ("entry_pdb_id", "entry_resolution")
APO_CANDIDATE_COLUMNS = (
    "target_system",
    "source_entry_id",
    "source_chain_asym_id",
    "source_chain_auth_id",
    "source_biounit_id",
    "source_chain_instance",
    "source_num_contacting_ions",
    "source_num_contacting_artifacts",
    "source_num_contacting_other_ligands",
    "source_resolution",
)


def _sql_path(path: str | Path) -> str:
    return Path(path).resolve().as_posix().replace("'", "''")


def _parquet_dataset(path: str | Path) -> ds.Dataset:
    source = Path(path)
    if source.is_dir():
        return ds.dataset(
            source,
            format="parquet",
            partitioning="hive",
            exclude_invalid_files=True,
        )
    return ds.dataset(source, format="parquet", exclude_invalid_files=True)


def _sql_parquet_source(path: str | Path) -> str:
    source = Path(path).resolve()
    if source.is_dir():
        source = source / "**/*.parquet"
    return source.as_posix().replace("'", "''")


def _check_parquet_columns(
    table: str | Path, columns: tuple[str, ...], name: str
) -> None:
    available = set(_parquet_dataset(table).schema.names)
    missing = sorted(set(columns).difference(available))
    if missing:
        raise ValueError(f"{name} is missing columns {missing}")


@dataclass(frozen=True)
class LinkedApoSelectionConfig:
    """Similarity thresholds and the maximum apo links per holo system."""

    max_per_system: int = 5
    min_pocket_fident: int = 95
    min_protein_fident_weighted_sum: int = 95
    min_protein_fident_qcov_weighted_sum: int = 80
    min_protein_lddt_weighted_sum: int = 20

    def __post_init__(self) -> None:
        if self.max_per_system < 1:
            raise ValueError("max_per_system must be positive")
        thresholds = {
            "min_pocket_fident": self.min_pocket_fident,
            "min_protein_fident_weighted_sum": (self.min_protein_fident_weighted_sum),
            "min_protein_fident_qcov_weighted_sum": (
                self.min_protein_fident_qcov_weighted_sum
            ),
            "min_protein_lddt_weighted_sum": self.min_protein_lddt_weighted_sum,
        }
        invalid = {
            name: value for name, value in thresholds.items() if not 0 <= value <= 100
        }
        if invalid:
            raise ValueError(f"linked-apo thresholds must be in [0, 100]: {invalid}")


def _read_columns(
    table: TableInput, columns: tuple[str, ...], name: str
) -> pd.DataFrame:
    if isinstance(table, pd.DataFrame):
        missing = sorted(set(columns).difference(table.columns))
        if missing:
            raise ValueError(f"{name} is missing columns {missing}")
        return table.loc[:, columns].copy()
    path = Path(table)
    available = set(_parquet_dataset(path).schema.names)
    missing = sorted(set(columns).difference(available))
    if missing:
        raise ValueError(f"{name} is missing columns {missing}")
    return pd.read_parquet(path, columns=list(columns))


def _read_score_columns(table: TableInput) -> pd.DataFrame:
    if isinstance(table, pd.DataFrame):
        return _read_columns(table, SCORE_COLUMNS, "protein score table")
    path = Path(table)
    available = set(_parquet_dataset(path).schema.names)
    missing = sorted(set(SCORE_COLUMNS).difference(available))
    if missing:
        raise ValueError(f"protein score table is missing columns {missing}")
    return pd.read_parquet(
        path,
        columns=list(SCORE_COLUMNS),
        filters=[("metric", "in", list(REQUIRED_SCORE_METRICS))],
    )


def _as_required_strings(frame: pd.DataFrame, columns: list[str], name: str) -> None:
    for column in columns:
        frame[column] = frame[column].astype("string")
        missing = frame[column].isna() | frame[column].str.strip().eq("")
        if missing.any():
            raise ValueError(f"{name} has empty {column} values")


def ligand_holo_chain_keys(annotation: TableInput) -> pd.DataFrame:
    """Return receptor-chain keys bound to proper ligands."""
    frame = _read_columns(
        annotation,
        HOLO_CHAIN_ANNOTATION_COLUMNS,
        "annotation table",
    )
    proper = frame["ligand_is_proper"].fillna(False).astype(bool)
    keys = frame.loc[proper].explode("ligand_protein_chains_asym_id", ignore_index=True)
    keys = keys.rename(columns={"ligand_protein_chains_asym_id": "chain_instance"})
    keys = keys[["entry_pdb_id", "chain_instance"]]
    encoded = keys["chain_instance"].astype("string")
    keys = keys.loc[encoded.notna() & encoded.str.strip().ne("")].copy()
    if keys.empty:
        return pd.DataFrame(columns=["entry_pdb_id", "chain_asym_id"])
    _as_required_strings(
        keys,
        ["entry_pdb_id", "chain_instance"],
        "annotation table proper-ligand chains",
    )
    keys["chain_asym_id"] = keys.pop("chain_instance").str.split(".", n=1).str[-1]
    return keys.drop_duplicates(ignore_index=True)


def build_apo_candidate_manifest(
    entry_chains: TableInput,
    *,
    biounit_chains: TableInput,
    entry_metadata: TableInput,
    annotation: TableInput,
) -> pd.DataFrame:
    """Build one scored target row per apo protein chain.

    A protein chain is apo when it is not ligand-like, is not bound to a
    proper ligand, and does not share an entity with a proper-ligand receptor
    chain in its entry. Protein-interface membership is independent of this
    ligand-relative definition. If a chain occurs in several biological
    assemblies, the assembly with the cleanest chain-local ligand environment
    is retained.
    """
    chains = _read_columns(entry_chains, ENTRY_CHAIN_COLUMNS, "entry chain table")
    membership = _read_columns(
        biounit_chains,
        BIOUNIT_CHAIN_COLUMNS,
        "biological-assembly membership table",
    )
    metadata = _read_columns(
        entry_metadata,
        ENTRY_METADATA_COLUMNS,
        "entry metadata table",
    )
    holo_keys = ligand_holo_chain_keys(annotation)
    _as_required_strings(
        chains,
        [
            "entry_pdb_id",
            "chain_asym_id",
            "chain_auth_id",
            "chain_receptor_type",
        ],
        "entry chain table",
    )
    chains["chain_entity_id"] = chains["chain_entity_id"].astype("string")
    duplicate_chains = chains.duplicated(["entry_pdb_id", "chain_asym_id"], keep=False)
    if duplicate_chains.any():
        examples = (
            chains.loc[duplicate_chains, ["entry_pdb_id", "chain_asym_id"]]
            .drop_duplicates()
            .head(10)
            .to_dict("records")
        )
        raise ValueError(f"entry chain table has duplicate chains: {examples}")

    holo_keys["chain_is_holo"] = True
    chains = chains.merge(
        holo_keys,
        on=["entry_pdb_id", "chain_asym_id"],
        how="left",
        validate="one_to_one",
    )
    holo = chains["chain_is_holo"].eq(True)
    usable_entity = chains["chain_entity_id"].notna() & chains[
        "chain_entity_id"
    ].str.strip().ne("")
    holo_entities = chains.loc[
        holo & usable_entity, ["entry_pdb_id", "chain_entity_id"]
    ].drop_duplicates()
    holo_entities["entity_is_holo"] = True
    ligand_like = chains["chain_is_ligand_like"].fillna(False).astype(bool)
    candidates = chains.loc[
        chains["chain_receptor_type"].str.lower().eq("protein") & ~holo & ~ligand_like
    ].merge(
        holo_entities,
        on=["entry_pdb_id", "chain_entity_id"],
        how="left",
        validate="many_to_one",
    )
    candidates = candidates.loc[candidates["entity_is_holo"].ne(True)].copy()
    if candidates.empty:
        return pd.DataFrame(columns=APO_CANDIDATE_COLUMNS)

    _as_required_strings(
        membership,
        [
            "entry_pdb_id",
            "biounit_id",
            "chain_instance",
            "chain_asym_id",
            "chain_role",
        ],
        "biological-assembly membership table",
    )
    membership_key = ["entry_pdb_id", "biounit_id", "chain_instance"]
    duplicate_membership = membership.duplicated(membership_key, keep=False)
    if duplicate_membership.any():
        examples = (
            membership.loc[duplicate_membership, membership_key]
            .drop_duplicates()
            .head(10)
            .to_dict("records")
        )
        raise ValueError(
            "biological-assembly membership has duplicate chain instances: "
            f"{examples}"
        )
    receptor_membership = (
        membership.loc[
            membership["chain_role"].str.lower().eq("receptor"),
            [
                "entry_pdb_id",
                "biounit_id",
                "chain_asym_id",
                "chain_instance",
                "chain_num_contacting_ions",
                "chain_num_contacting_artifacts",
                "chain_num_contacting_other_ligands",
            ],
        ]
        .sort_values("chain_instance")
        .rename(
            columns={
                "chain_num_contacting_ions": "source_num_contacting_ions",
                "chain_num_contacting_artifacts": ("source_num_contacting_artifacts"),
                "chain_num_contacting_other_ligands": (
                    "source_num_contacting_other_ligands"
                ),
            }
        )
    )
    candidates = candidates.merge(
        receptor_membership,
        on=["entry_pdb_id", "chain_asym_id"],
        how="inner",
        validate="one_to_many",
    )
    contact_columns = [
        "source_num_contacting_ions",
        "source_num_contacting_artifacts",
        "source_num_contacting_other_ligands",
    ]
    for column in contact_columns:
        candidates[column] = pd.to_numeric(candidates[column], errors="coerce")
        invalid = ~(
            np.isfinite(candidates[column])
            & candidates[column].ge(0)
            & candidates[column].mod(1).eq(0)
        )
        if invalid.any():
            raise ValueError(f"biological-assembly membership has invalid {column}")
        candidates[column] = candidates[column].astype("int64")

    _as_required_strings(metadata, ["entry_pdb_id"], "entry metadata table")
    duplicate_metadata = metadata["entry_pdb_id"].duplicated(keep=False)
    if duplicate_metadata.any():
        duplicate_ids = sorted(
            metadata.loc[duplicate_metadata, "entry_pdb_id"].unique()
        )
        raise ValueError(
            f"entry metadata table has duplicate entries: {duplicate_ids[:10]}"
        )
    metadata["entry_resolution"] = pd.to_numeric(
        metadata["entry_resolution"], errors="coerce"
    )
    invalid_resolution = metadata["entry_resolution"].notna() & ~(
        np.isfinite(metadata["entry_resolution"]) & metadata["entry_resolution"].gt(0)
    )
    if invalid_resolution.any():
        invalid = metadata.loc[
            invalid_resolution, ["entry_pdb_id", "entry_resolution"]
        ].to_dict("records")
        raise ValueError(f"entry metadata has invalid resolutions: {invalid[:10]}")
    candidates = candidates.merge(
        metadata,
        on="entry_pdb_id",
        how="left",
        validate="many_to_one",
    )
    candidates["resolution_missing"] = candidates["entry_resolution"].isna()
    candidates["ligand_contact_class"] = np.select(
        [
            candidates["source_num_contacting_other_ligands"].gt(0),
            candidates["source_num_contacting_artifacts"].gt(0),
            candidates["source_num_contacting_ions"].gt(0),
        ],
        [3, 2, 1],
        default=0,
    )
    candidates["num_contacting_ligands"] = candidates[contact_columns].sum(axis=1)
    candidates = candidates.sort_values(
        [
            "entry_pdb_id",
            "chain_asym_id",
            "ligand_contact_class",
            "num_contacting_ligands",
            "resolution_missing",
            "entry_resolution",
            "biounit_id",
            "chain_instance",
        ],
        ignore_index=True,
    ).drop_duplicates(["entry_pdb_id", "chain_asym_id"], keep="first")
    candidates["target_system"] = (
        candidates["entry_pdb_id"] + "_" + candidates["chain_asym_id"]
    )
    candidates = candidates.rename(
        columns={
            "entry_pdb_id": "source_entry_id",
            "chain_asym_id": "source_chain_asym_id",
            "chain_auth_id": "source_chain_auth_id",
            "biounit_id": "source_biounit_id",
            "chain_instance": "source_chain_instance",
            "entry_resolution": "source_resolution",
        }
    )
    return candidates.loc[:, APO_CANDIDATE_COLUMNS].reset_index(drop=True)


def _validate_apo_candidates(candidates: pd.DataFrame) -> pd.DataFrame:
    candidates = candidates.copy()
    _as_required_strings(
        candidates,
        [
            "target_system",
            "source_entry_id",
            "source_chain_asym_id",
            "source_chain_auth_id",
            "source_biounit_id",
            "source_chain_instance",
        ],
        "apo candidate manifest",
    )
    duplicate_ids = sorted(
        candidates.loc[
            candidates["target_system"].duplicated(keep=False), "target_system"
        ].unique()
    )
    if duplicate_ids:
        raise ValueError(
            "apo candidate manifest has duplicate target_system values: "
            f"{duplicate_ids[:10]}"
        )
    for column in [
        "source_num_contacting_ions",
        "source_num_contacting_artifacts",
        "source_num_contacting_other_ligands",
    ]:
        candidates[column] = pd.to_numeric(candidates[column], errors="coerce")
        invalid_counts = ~(
            np.isfinite(candidates[column])
            & candidates[column].ge(0)
            & candidates[column].mod(1).eq(0)
        )
        if invalid_counts.any():
            raise ValueError(
                "apo candidate contact counts must be non-negative integers"
            )
        candidates[column] = candidates[column].astype("int64")
    candidates["source_resolution"] = pd.to_numeric(
        candidates["source_resolution"], errors="coerce"
    )
    invalid_resolution = candidates["source_resolution"].notna() & ~(
        np.isfinite(candidates["source_resolution"])
        & candidates["source_resolution"].gt(0)
    )
    if invalid_resolution.any():
        raise ValueError("apo candidate resolutions must be positive when present")
    return candidates


def _proper_ligands(annotation: pd.DataFrame) -> pd.DataFrame:
    proper = annotation.loc[
        annotation["ligand_is_proper"].fillna(False).astype(bool),
        ["system_id", "ligand_id"],
    ].copy()
    _as_required_strings(proper, ["system_id", "ligand_id"], "proper ligand rows")
    return proper.drop_duplicates(ignore_index=True)


def _prepare_scores(scores: pd.DataFrame) -> pd.DataFrame:
    scores = scores.copy()
    _as_required_strings(
        scores,
        ["query_system", "query_ligand_id", "target_system", "metric"],
        "protein score table",
    )
    scores = scores.loc[scores["metric"].isin(REQUIRED_SCORE_METRICS)].copy()
    scores["similarity"] = pd.to_numeric(scores["similarity"], errors="coerce")
    invalid_similarity = ~(
        np.isfinite(scores["similarity"])
        & scores["similarity"].between(0, 100, inclusive="both")
    )
    if invalid_similarity.any():
        raise ValueError(
            "required protein similarities must be finite values in [0, 100]"
        )
    key = ["query_system", "query_ligand_id", "target_system", "metric"]
    duplicate = scores.duplicated(key, keep=False)
    if duplicate.any():
        examples = (
            scores.loc[duplicate, key].drop_duplicates().head(10).to_dict("records")
        )
        raise ValueError(
            f"protein score table has duplicate required metrics: {examples}"
        )
    return scores


def select_linked_apo_structures(
    protein_scores: TableInput,
    *,
    annotation: TableInput,
    candidates: TableInput,
    config: LinkedApoSelectionConfig | None = None,
) -> pd.DataFrame:
    """Associate PLINDER holo systems with matching deposited apo chains.

    A candidate must satisfy every required metric for every proper ligand
    pocket in the holo system. Candidates from the holo entry itself are
    excluded. Passing candidates prefer no chain-local ligand contacts, then
    ion-only contacts, artifact contacts, and other ligand contacts before
    experimental resolution and similarity.
    """
    config = config or LinkedApoSelectionConfig()
    score_frame = _prepare_scores(_read_score_columns(protein_scores))
    candidate_frame = _validate_apo_candidates(
        _read_columns(candidates, APO_CANDIDATE_COLUMNS, "apo candidate manifest")
    )
    proper = _proper_ligands(
        _read_columns(annotation, ANNOTATION_COLUMNS, "annotation table")
    )
    if score_frame.empty or candidate_frame.empty or proper.empty:
        return pd.DataFrame(columns=STRUCTURE_LINK_SCHEMA.names)

    score_systems = set(score_frame["query_system"])
    known_systems = set(proper["system_id"])
    unknown_systems = sorted(score_systems.difference(known_systems))
    if unknown_systems:
        raise ValueError(
            "protein scores reference systems with no proper ligand rows: "
            f"{unknown_systems[:10]}"
        )

    scored = score_frame.merge(
        candidate_frame,
        on="target_system",
        how="inner",
        validate="many_to_one",
    ).merge(
        proper,
        left_on=["query_system", "query_ligand_id"],
        right_on=["system_id", "ligand_id"],
        how="inner",
        validate="many_to_one",
    )
    scored = scored.loc[
        scored["query_system"].str.split("__", n=1).str[0] != scored["source_entry_id"]
    ]
    if scored.empty:
        return pd.DataFrame(columns=STRUCTURE_LINK_SCHEMA.names)

    ligand_metrics = scored.pivot(
        index=["query_system", "query_ligand_id", "target_system"],
        columns="metric",
        values="similarity",
    ).reset_index()
    for metric in REQUIRED_SCORE_METRICS:
        if metric not in ligand_metrics:
            ligand_metrics[metric] = np.nan
    ligand_metrics = ligand_metrics.dropna(subset=list(REQUIRED_SCORE_METRICS))

    expected_ligands = (
        proper.groupby("system_id", observed=True)["ligand_id"]
        .nunique()
        .rename("expected_ligand_pockets")
    )
    summaries = (
        ligand_metrics.groupby(
            ["query_system", "target_system"], observed=True, as_index=False
        )
        .agg(
            num_ligand_pockets=("query_ligand_id", "nunique"),
            min_pocket_fident=("pocket_fident", "min"),
            mean_pocket_fident=("pocket_fident", "mean"),
            min_protein_fident_weighted_sum=(
                "protein_fident_weighted_sum",
                "min",
            ),
            min_protein_fident_qcov_weighted_sum=(
                "protein_fident_qcov_weighted_sum",
                "min",
            ),
            min_protein_lddt_weighted_sum=(
                "protein_lddt_weighted_sum",
                "min",
            ),
        )
        .merge(
            expected_ligands,
            left_on="query_system",
            right_index=True,
            how="left",
            validate="many_to_one",
        )
    )
    summaries = summaries.loc[
        (summaries["num_ligand_pockets"] == summaries["expected_ligand_pockets"])
        & (summaries["min_pocket_fident"] >= config.min_pocket_fident)
        & (
            summaries["min_protein_fident_weighted_sum"]
            >= config.min_protein_fident_weighted_sum
        )
        & (
            summaries["min_protein_fident_qcov_weighted_sum"]
            >= config.min_protein_fident_qcov_weighted_sum
        )
        & (
            summaries["min_protein_lddt_weighted_sum"]
            >= config.min_protein_lddt_weighted_sum
        )
    ].copy()
    if summaries.empty:
        return pd.DataFrame(columns=STRUCTURE_LINK_SCHEMA.names)

    links = summaries.merge(
        candidate_frame,
        on="target_system",
        how="inner",
        validate="many_to_one",
    ).rename(
        columns={
            "query_system": "reference_system_id",
            "target_system": "linked_structure_id",
        }
    )
    links["resolution_missing"] = links["source_resolution"].isna()
    links["ligand_contact_class"] = np.select(
        [
            links["source_num_contacting_other_ligands"].gt(0),
            links["source_num_contacting_artifacts"].gt(0),
            links["source_num_contacting_ions"].gt(0),
        ],
        [3, 2, 1],
        default=0,
    )
    links["num_contacting_ligands"] = links[
        [
            "source_num_contacting_ions",
            "source_num_contacting_artifacts",
            "source_num_contacting_other_ligands",
        ]
    ].sum(axis=1)
    links = links.sort_values(
        [
            "reference_system_id",
            "ligand_contact_class",
            "num_contacting_ligands",
            "resolution_missing",
            "source_resolution",
            "min_pocket_fident",
            "mean_pocket_fident",
            "min_protein_fident_qcov_weighted_sum",
            "source_entry_id",
            "source_chain_asym_id",
        ],
        ascending=[
            True,
            True,
            True,
            True,
            True,
            False,
            False,
            False,
            True,
            True,
        ],
        ignore_index=True,
    )
    links["rank"] = links.groupby("reference_system_id", observed=True).cumcount() + 1
    return links.loc[
        links["rank"] <= config.max_per_system,
        STRUCTURE_LINK_SCHEMA.names,
    ].reset_index(drop=True)


def write_linked_apo_structure_table(
    protein_scores: TableInput,
    *,
    annotation: TableInput,
    candidates: TableInput,
    output_path: str | Path,
    config: LinkedApoSelectionConfig | None = None,
    scratch_dir: str | Path | None = None,
    threads: int = 1,
    memory_limit: str = "7GB",
) -> Path:
    """Select linked apo structures and write the compact release table."""
    output_path = Path(output_path)
    if all(
        not isinstance(table, pd.DataFrame)
        for table in (protein_scores, annotation, candidates)
    ):
        return _write_linked_apo_structure_table_from_parquet(
            protein_scores=protein_scores,
            annotation=annotation,
            candidates=candidates,
            output_path=output_path,
            config=config or LinkedApoSelectionConfig(),
            scratch_dir=scratch_dir,
            threads=threads,
            memory_limit=memory_limit,
        )
    links = select_linked_apo_structures(
        protein_scores,
        annotation=annotation,
        candidates=candidates,
        config=config,
    )
    output_path.parent.mkdir(exist_ok=True, parents=True)
    table = pa.Table.from_pandas(
        links,
        schema=STRUCTURE_LINK_SCHEMA,
        preserve_index=False,
        safe=True,
    )
    pq.write_table(table, output_path, compression="zstd")
    return output_path


def _write_linked_apo_structure_table_from_parquet(
    *,
    protein_scores: str | Path,
    annotation: str | Path,
    candidates: str | Path,
    output_path: Path,
    config: LinkedApoSelectionConfig,
    scratch_dir: str | Path | None,
    threads: int,
    memory_limit: str,
) -> Path:
    """Aggregate a release-scale score table without loading it into Pandas."""
    import duckdb

    if threads < 1:
        raise ValueError("threads must be positive")
    for table, columns, name in (
        (protein_scores, SCORE_COLUMNS, "protein score table"),
        (annotation, ANNOTATION_COLUMNS, "annotation table"),
        (candidates, APO_CANDIDATE_COLUMNS, "apo candidate manifest"),
    ):
        _check_parquet_columns(table, columns, name)

    scratch = Path(scratch_dir or output_path.parent / ".linked-apo-scratch")
    scratch.mkdir(exist_ok=True, parents=True)
    output_path.parent.mkdir(exist_ok=True, parents=True)
    temporary = output_path.with_suffix(output_path.suffix + ".tmp")
    temporary.unlink(missing_ok=True)
    metrics_sql = ", ".join(f"'{metric}'" for metric in REQUIRED_SCORE_METRICS)

    connection = duckdb.connect()
    try:
        connection.sql(f"SET threads={threads}")
        connection.sql(f"SET temp_directory='{_sql_path(scratch)}'")
        connection.sql(f"SET memory_limit='{memory_limit}'")
        connection.sql("SET preserve_insertion_order=false")
        duplicate_scores = connection.sql(
            f"""
            SELECT
                CAST(query_system AS VARCHAR) AS query_system,
                CAST(query_ligand_id AS VARCHAR) AS query_ligand_id,
                CAST(target_system AS VARCHAR) AS target_system,
                CAST(metric AS VARCHAR) AS metric
            FROM read_parquet('{_sql_parquet_source(protein_scores)}')
            WHERE metric IN ({metrics_sql})
            GROUP BY query_system, query_ligand_id, target_system, metric
            HAVING COUNT(*) > 1
            LIMIT 10
            """
        ).df()
        if not duplicate_scores.empty:
            raise ValueError(
                "protein score table has duplicate required metrics: "
                f"{duplicate_scores.to_dict('records')}"
            )
        if _parquet_dataset(candidates).count_rows() == 0:
            pq.write_table(
                pa.Table.from_pylist([], schema=STRUCTURE_LINK_SCHEMA),
                temporary,
                compression="zstd",
            )
            temporary.replace(output_path)
            return output_path
        connection.sql(
            f"""
            COPY (
                WITH proper AS (
                    SELECT DISTINCT
                        CAST(system_id AS VARCHAR) AS system_id,
                        CAST(ligand_id AS VARCHAR) AS ligand_id
                    FROM read_parquet('{_sql_parquet_source(annotation)}')
                    WHERE COALESCE(
                        TRY_CAST(ligand_is_proper AS BOOLEAN), FALSE
                    )
                ),
                expected AS (
                    SELECT
                        system_id,
                        COUNT(DISTINCT ligand_id) AS expected_ligand_pockets
                    FROM proper
                    GROUP BY system_id
                ),
                score_base AS (
                    SELECT
                        CAST(query_system AS VARCHAR) AS query_system,
                        CAST(query_ligand_id AS VARCHAR) AS query_ligand_id,
                        CAST(target_system AS VARCHAR) AS target_system,
                        CAST(metric AS VARCHAR) AS metric,
                        TRY_CAST(similarity AS DOUBLE) AS similarity
                    FROM read_parquet('{_sql_parquet_source(protein_scores)}')
                    WHERE metric IN ({metrics_sql})
                ),
                scored AS (
                    SELECT
                        scores.query_system,
                        scores.query_ligand_id,
                        scores.target_system,
                        scores.metric,
                        scores.similarity
                    FROM score_base AS scores
                    INNER JOIN read_parquet(
                        '{_sql_parquet_source(candidates)}'
                    ) AS candidate USING (target_system)
                    INNER JOIN proper
                        ON scores.query_system = proper.system_id
                        AND scores.query_ligand_id = proper.ligand_id
                    WHERE split_part(scores.query_system, '__', 1)
                        != candidate.source_entry_id
                        AND scores.similarity BETWEEN 0 AND 100
                ),
                per_ligand AS (
                    SELECT
                        query_system,
                        query_ligand_id,
                        target_system,
                        MAX(similarity) FILTER (
                            WHERE metric = 'pocket_fident'
                        ) AS pocket_fident,
                        MAX(similarity) FILTER (
                            WHERE metric = 'protein_fident_weighted_sum'
                        ) AS protein_fident_weighted_sum,
                        MAX(similarity) FILTER (
                            WHERE metric = 'protein_fident_qcov_weighted_sum'
                        ) AS protein_fident_qcov_weighted_sum,
                        MAX(similarity) FILTER (
                            WHERE metric = 'protein_lddt_weighted_sum'
                        ) AS protein_lddt_weighted_sum
                    FROM scored
                    GROUP BY query_system, query_ligand_id, target_system
                    HAVING COUNT(*) = 4 AND COUNT(DISTINCT metric) = 4
                ),
                summaries AS (
                    SELECT
                        query_system,
                        target_system,
                        COUNT(DISTINCT query_ligand_id) AS num_ligand_pockets,
                        MIN(pocket_fident) AS min_pocket_fident,
                        AVG(pocket_fident) AS mean_pocket_fident,
                        MIN(protein_fident_weighted_sum)
                            AS min_protein_fident_weighted_sum,
                        MIN(protein_fident_qcov_weighted_sum)
                            AS min_protein_fident_qcov_weighted_sum,
                        MIN(protein_lddt_weighted_sum)
                            AS min_protein_lddt_weighted_sum
                    FROM per_ligand
                    GROUP BY query_system, target_system
                ),
                eligible AS (
                    SELECT summaries.*, candidate.* EXCLUDE (target_system)
                    FROM summaries
                    INNER JOIN expected
                        ON summaries.query_system = expected.system_id
                    INNER JOIN read_parquet(
                        '{_sql_parquet_source(candidates)}'
                    ) AS candidate USING (target_system)
                    WHERE summaries.num_ligand_pockets
                            = expected.expected_ligand_pockets
                        AND summaries.min_pocket_fident
                            >= {config.min_pocket_fident}
                        AND summaries.min_protein_fident_weighted_sum
                            >= {config.min_protein_fident_weighted_sum}
                        AND summaries.min_protein_fident_qcov_weighted_sum
                            >= {config.min_protein_fident_qcov_weighted_sum}
                        AND summaries.min_protein_lddt_weighted_sum
                            >= {config.min_protein_lddt_weighted_sum}
                ),
                ranked AS (
                    SELECT
                        *,
                        ROW_NUMBER() OVER (
                            PARTITION BY query_system
                            ORDER BY
                                CASE
                                    WHEN source_num_contacting_other_ligands > 0
                                        THEN 3
                                    WHEN source_num_contacting_artifacts > 0
                                        THEN 2
                                    WHEN source_num_contacting_ions > 0 THEN 1
                                    ELSE 0
                                END,
                                source_num_contacting_ions
                                    + source_num_contacting_artifacts
                                    + source_num_contacting_other_ligands,
                                source_resolution IS NULL,
                                source_resolution,
                                min_pocket_fident DESC,
                                mean_pocket_fident DESC,
                                min_protein_fident_qcov_weighted_sum DESC,
                                source_entry_id,
                                source_chain_asym_id
                        ) AS link_rank
                    FROM eligible
                )
                SELECT
                    CAST(query_system AS VARCHAR) AS reference_system_id,
                    CAST(target_system AS VARCHAR) AS linked_structure_id,
                    CAST(source_entry_id AS VARCHAR) AS source_entry_id,
                    CAST(source_chain_asym_id AS VARCHAR)
                        AS source_chain_asym_id,
                    CAST(source_chain_auth_id AS VARCHAR)
                        AS source_chain_auth_id,
                    CAST(source_biounit_id AS VARCHAR) AS source_biounit_id,
                    CAST(source_chain_instance AS VARCHAR)
                        AS source_chain_instance,
                    CAST(source_num_contacting_ions AS SMALLINT)
                        AS source_num_contacting_ions,
                    CAST(source_num_contacting_artifacts AS SMALLINT)
                        AS source_num_contacting_artifacts,
                    CAST(source_num_contacting_other_ligands AS SMALLINT)
                        AS source_num_contacting_other_ligands,
                    CAST(source_resolution AS FLOAT) AS source_resolution,
                    CAST(link_rank AS SMALLINT) AS rank,
                    CAST(num_ligand_pockets AS SMALLINT) AS num_ligand_pockets,
                    CAST(min_pocket_fident AS TINYINT) AS min_pocket_fident,
                    CAST(mean_pocket_fident AS FLOAT) AS mean_pocket_fident,
                    CAST(min_protein_fident_weighted_sum AS TINYINT)
                        AS min_protein_fident_weighted_sum,
                    CAST(min_protein_fident_qcov_weighted_sum AS TINYINT)
                        AS min_protein_fident_qcov_weighted_sum,
                    CAST(min_protein_lddt_weighted_sum AS TINYINT)
                        AS min_protein_lddt_weighted_sum
                FROM ranked
                WHERE link_rank <= {config.max_per_system}
                ORDER BY reference_system_id, rank
            ) TO '{_sql_path(temporary)}' (
                FORMAT PARQUET,
                COMPRESSION ZSTD,
                ROW_GROUP_SIZE 500000
            )
            """
        )
        temporary.replace(output_path)
    finally:
        connection.close()
        temporary.unlink(missing_ok=True)
    return output_path
