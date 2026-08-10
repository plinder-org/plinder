# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import TypeAlias

import numpy as np
import pandas as pd
import pyarrow as pa
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
ENTRY_CHAIN_COLUMNS = (
    "entry_pdb_id",
    "chain_asym_id",
    "chain_auth_id",
    "chain_entity_id",
    "chain_receptor_type",
    "chain_is_holo",
)
BIOUNIT_CHAIN_COLUMNS = (
    "entry_pdb_id",
    "biounit_id",
    "chain_instance",
    "chain_asym_id",
    "chain_role",
)
ENTRY_METADATA_COLUMNS = ("entry_pdb_id", "entry_resolution")
APO_CANDIDATE_COLUMNS = (
    "target_system",
    "source_entry_id",
    "source_chain_asym_id",
    "source_chain_auth_id",
    "source_biounit_id",
    "source_chain_instance",
    "source_num_ligand_chains",
    "source_resolution",
)


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
            "min_protein_fident_weighted_sum": (
                self.min_protein_fident_weighted_sum
            ),
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


def _read_columns(table: TableInput, columns: tuple[str, ...], name: str) -> pd.DataFrame:
    if isinstance(table, pd.DataFrame):
        missing = sorted(set(columns).difference(table.columns))
        if missing:
            raise ValueError(f"{name} is missing columns {missing}")
        return table.loc[:, columns].copy()
    path = Path(table)
    available = set(pq.read_schema(path).names)
    missing = sorted(set(columns).difference(available))
    if missing:
        raise ValueError(f"{name} is missing columns {missing}")
    return pd.read_parquet(path, columns=list(columns))


def _read_score_columns(table: TableInput) -> pd.DataFrame:
    if isinstance(table, pd.DataFrame):
        return _read_columns(table, SCORE_COLUMNS, "protein score table")
    path = Path(table)
    available = set(pq.read_schema(path).names)
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


def build_apo_candidate_manifest(
    entry_chains: TableInput,
    *,
    biounit_chains: TableInput,
    entry_metadata: TableInput,
) -> pd.DataFrame:
    """Build one scored target row per apo protein chain.

    The apo definition matches the alignment database: a protein chain must
    not be holo and must not share an entity with a holo chain in its entry.
    If the chain occurs in several biological assemblies, the assembly with
    the fewest non-water ligand chains is retained.
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
    duplicate_chains = chains.duplicated(
        ["entry_pdb_id", "chain_asym_id"], keep=False
    )
    if duplicate_chains.any():
        examples = (
            chains.loc[duplicate_chains, ["entry_pdb_id", "chain_asym_id"]]
            .drop_duplicates()
            .head(10)
            .to_dict("records")
        )
        raise ValueError(f"entry chain table has duplicate chains: {examples}")

    holo = chains["chain_is_holo"].fillna(False).astype(bool)
    usable_entity = chains["chain_entity_id"].notna() & chains[
        "chain_entity_id"
    ].str.strip().ne("")
    holo_entities = chains.loc[
        holo & usable_entity, ["entry_pdb_id", "chain_entity_id"]
    ].drop_duplicates()
    holo_entities["entity_is_holo"] = True
    candidates = chains.loc[
        chains["chain_receptor_type"].str.lower().eq("protein") & ~holo
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
    ligand_counts = (
        membership.loc[membership["chain_role"].str.lower().eq("ligand")]
        .groupby(["entry_pdb_id", "biounit_id"], observed=True)["chain_instance"]
        .nunique()
        .rename("source_num_ligand_chains")
        .reset_index()
    )
    receptor_membership = (
        membership.loc[
            membership["chain_role"].str.lower().eq("receptor"),
            ["entry_pdb_id", "biounit_id", "chain_asym_id", "chain_instance"],
        ]
        .sort_values("chain_instance")
        .drop_duplicates(["entry_pdb_id", "biounit_id", "chain_asym_id"])
    )
    candidates = candidates.merge(
        receptor_membership,
        on=["entry_pdb_id", "chain_asym_id"],
        how="left",
        validate="one_to_many",
    )
    missing_membership = candidates["biounit_id"].isna()
    if missing_membership.any():
        missing = candidates.loc[
            missing_membership, ["entry_pdb_id", "chain_asym_id"]
        ].to_dict("records")
        raise ValueError(f"apo chains have no biological-assembly membership: {missing[:10]}")
    candidates = candidates.merge(
        ligand_counts,
        on=["entry_pdb_id", "biounit_id"],
        how="left",
        validate="many_to_one",
    )
    candidates["source_num_ligand_chains"] = (
        candidates["source_num_ligand_chains"].fillna(0).astype("int64")
    )

    _as_required_strings(metadata, ["entry_pdb_id"], "entry metadata table")
    duplicate_metadata = metadata["entry_pdb_id"].duplicated(keep=False)
    if duplicate_metadata.any():
        duplicate_ids = sorted(metadata.loc[duplicate_metadata, "entry_pdb_id"].unique())
        raise ValueError(f"entry metadata table has duplicate entries: {duplicate_ids[:10]}")
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
    candidates = candidates.sort_values(
        [
            "entry_pdb_id",
            "chain_asym_id",
            "source_num_ligand_chains",
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
    candidates["source_num_ligand_chains"] = pd.to_numeric(
        candidates["source_num_ligand_chains"], errors="coerce"
    )
    invalid_counts = ~(
        np.isfinite(candidates["source_num_ligand_chains"])
        & candidates["source_num_ligand_chains"].ge(0)
        & candidates["source_num_ligand_chains"].mod(1).eq(0)
    )
    if invalid_counts.any():
        raise ValueError("apo candidate ligand-chain counts must be non-negative integers")
    candidates["source_num_ligand_chains"] = candidates[
        "source_num_ligand_chains"
    ].astype("int64")
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
        raise ValueError("required protein similarities must be finite values in [0, 100]")
    key = ["query_system", "query_ligand_id", "target_system", "metric"]
    duplicate = scores.duplicated(key, keep=False)
    if duplicate.any():
        examples = scores.loc[duplicate, key].drop_duplicates().head(10).to_dict("records")
        raise ValueError(f"protein score table has duplicate required metrics: {examples}")
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
    excluded. Passing candidates are ranked by assembly ligand count first,
    then by experimental resolution and similarity.
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
        scored["query_system"].str.split("__", n=1).str[0]
        != scored["source_entry_id"]
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
    links = links.sort_values(
        [
            "reference_system_id",
            "source_num_ligand_chains",
            "resolution_missing",
            "source_resolution",
            "min_pocket_fident",
            "mean_pocket_fident",
            "min_protein_fident_qcov_weighted_sum",
            "source_entry_id",
            "source_chain_asym_id",
        ],
        ascending=[True, True, True, True, False, False, False, True, True],
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
) -> Path:
    """Select linked apo structures and write the compact release table."""
    links = select_linked_apo_structures(
        protein_scores,
        annotation=annotation,
        candidates=candidates,
        config=config,
    )
    output_path = Path(output_path)
    output_path.parent.mkdir(exist_ok=True, parents=True)
    table = pa.Table.from_pandas(
        links,
        schema=STRUCTURE_LINK_SCHEMA,
        preserve_index=False,
        safe=True,
    )
    pq.write_table(table, output_path, compression="zstd")
    return output_path
