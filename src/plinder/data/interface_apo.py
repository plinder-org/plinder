# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
"""Pair each protein-interface side with deposited non-interacting chains."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, TypeAlias

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.dataset as ds
import pyarrow.parquet as pq

from plinder.core.utils.schemas import INTERFACE_APO_LINK_SCHEMA

TableInput: TypeAlias = pd.DataFrame | str | Path

INTERFACE_COLUMNS = (
    "entry_pdb_id",
    "system_id",
    "interface_chain_1",
    "interface_chain_2",
)
ENTRY_CHAIN_COLUMNS = (
    "entry_pdb_id",
    "chain_asym_id",
    "chain_auth_id",
    "chain_receptor_type",
    "chain_is_ligand_like",
)
BIOUNIT_CHAIN_COLUMNS = (
    "entry_pdb_id",
    "biounit_id",
    "chain_instance",
    "chain_asym_id",
    "chain_role",
    "chain_num_contacting_proteins",
    "chain_num_contacting_ions",
    "chain_num_contacting_artifacts",
    "chain_num_contacting_other_ligands",
)
ENTRY_METADATA_COLUMNS = ("entry_pdb_id", "entry_resolution")
QUERY_COLUMNS = (
    "query_entry",
    "reference_system_id",
    "reference_side",
    "reference_chain_instance",
    "reference_chain_asym_id",
)
CANDIDATE_COLUMNS = (
    "target_system",
    "source_entry_id",
    "source_chain_asym_id",
    "source_chain_auth_id",
    "source_biounit_id",
    "source_chain_instance",
    "source_num_contacting_proteins",
    "source_num_contacting_ions",
    "source_num_contacting_artifacts",
    "source_num_contacting_other_ligands",
    "source_resolution",
)
ALIGNMENT_COLUMNS = (
    "query_entry",
    "target_entry",
    "query_chain_mapped",
    "target_chain_mapped",
    "qcov",
    "tcov",
    "fident",
)


@dataclass(frozen=True)
class InterfaceApoSelectionConfig:
    """Similarity thresholds and maximum links for each interface side."""

    max_per_side: int = 5
    min_sequence_identity: float = 0.95
    min_query_coverage: float = 0.80
    min_target_coverage: float = 0.80

    def __post_init__(self) -> None:
        if self.max_per_side < 1:
            raise ValueError("max_per_side must be positive")
        thresholds = {
            "min_sequence_identity": self.min_sequence_identity,
            "min_query_coverage": self.min_query_coverage,
            "min_target_coverage": self.min_target_coverage,
        }
        invalid = {
            name: value for name, value in thresholds.items() if not 0 <= value <= 1
        }
        if invalid:
            raise ValueError(f"interface apo thresholds must be in [0, 1]: {invalid}")


def _dataset(path: str | Path) -> Any:
    source = Path(path)
    if source.is_dir():
        return ds.dataset(
            source,
            format="parquet",
            partitioning="hive",
            exclude_invalid_files=True,
        )
    return ds.dataset(source, format="parquet", exclude_invalid_files=True)


def _read_columns(
    table: TableInput, columns: tuple[str, ...], name: str
) -> pd.DataFrame:
    if isinstance(table, pd.DataFrame):
        available = set(table.columns)
        frame = table
    else:
        available = set(_dataset(table).schema.names)
    missing = sorted(set(columns).difference(available))
    if missing:
        raise ValueError(f"{name} is missing columns {missing}")
    if not isinstance(table, pd.DataFrame):
        frame = pd.read_parquet(table, columns=list(columns))
    return frame.loc[:, columns].copy()


def _required_strings(frame: pd.DataFrame, columns: tuple[str, ...], name: str) -> None:
    for column in columns:
        frame[column] = frame[column].astype("string")
        if (frame[column].isna() | frame[column].str.strip().eq("")).any():
            raise ValueError(f"{name} has empty {column} values")


def build_interface_apo_query_manifest(interfaces: TableInput) -> pd.DataFrame:
    """Return one query row for each side of every protein interface."""
    frame = _read_columns(interfaces, INTERFACE_COLUMNS, "interface table")
    _required_strings(frame, INTERFACE_COLUMNS, "interface table")
    duplicate = frame["system_id"].duplicated(keep=False)
    if duplicate.any():
        examples = sorted(frame.loc[duplicate, "system_id"].unique())[:10]
        raise ValueError(f"interface table has duplicate systems: {examples}")

    sides = []
    for side in (1, 2):
        side_frame = frame[
            ["entry_pdb_id", "system_id", f"interface_chain_{side}"]
        ].rename(
            columns={
                "entry_pdb_id": "query_entry",
                "system_id": "reference_system_id",
                f"interface_chain_{side}": "reference_chain_instance",
            }
        )
        side_frame["reference_side"] = side
        sides.append(side_frame)
    queries = pd.concat(sides, ignore_index=True)
    queries["reference_chain_asym_id"] = (
        queries["reference_chain_instance"].str.split(".", n=1).str[-1]
    )
    return queries.loc[:, QUERY_COLUMNS].sort_values(
        ["query_entry", "reference_system_id", "reference_side"],
        ignore_index=True,
    )


def _contact_counts(
    frame: pd.DataFrame, columns: list[str], *, allow_missing: bool = False
) -> None:
    for column in columns:
        frame[column] = pd.to_numeric(frame[column], errors="coerce")
        valid = (
            np.isfinite(frame[column])
            & frame[column].ge(0)
            & frame[column].mod(1).eq(0)
        )
        if allow_missing:
            valid |= frame[column].isna()
        if not valid.fillna(False).all():
            raise ValueError(f"biological-assembly membership has invalid {column}")
        frame[column] = frame[column].astype("Int64" if allow_missing else "int64")


def build_interface_apo_candidate_manifest(
    entry_chains: TableInput,
    *,
    biounit_chains: TableInput,
    entry_metadata: TableInput,
) -> pd.DataFrame:
    """Choose a protein-contact-free assembly instance for each candidate chain."""
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
    _required_strings(
        chains,
        ("entry_pdb_id", "chain_asym_id", "chain_auth_id", "chain_receptor_type"),
        "entry chain table",
    )
    if chains.duplicated(["entry_pdb_id", "chain_asym_id"]).any():
        raise ValueError("entry chain table has duplicate entry/asym chain keys")
    chains = chains.loc[
        chains["chain_receptor_type"].str.lower().eq("protein")
        & ~chains["chain_is_ligand_like"].fillna(False).astype(bool)
    ].copy()

    _required_strings(
        membership,
        (
            "entry_pdb_id",
            "biounit_id",
            "chain_instance",
            "chain_asym_id",
            "chain_role",
        ),
        "biological-assembly membership table",
    )
    if membership.duplicated(["entry_pdb_id", "biounit_id", "chain_instance"]).any():
        raise ValueError("biological-assembly membership has duplicate chain instances")
    membership = membership.loc[
        membership["chain_role"].str.lower().eq("receptor")
    ].copy()
    candidates = chains.merge(
        membership,
        on=["entry_pdb_id", "chain_asym_id"],
        how="inner",
        validate="one_to_many",
    )
    ligand_contact_columns = [
        "chain_num_contacting_ions",
        "chain_num_contacting_artifacts",
        "chain_num_contacting_other_ligands",
    ]
    _contact_counts(candidates, ["chain_num_contacting_proteins"])
    _contact_counts(candidates, ligand_contact_columns, allow_missing=True)
    candidates = candidates.loc[candidates["chain_num_contacting_proteins"].eq(0)]
    if candidates.empty:
        return pd.DataFrame(columns=CANDIDATE_COLUMNS)

    _required_strings(metadata, ("entry_pdb_id",), "entry metadata table")
    if metadata["entry_pdb_id"].duplicated().any():
        raise ValueError("entry metadata table has duplicate entries")
    metadata["entry_resolution"] = pd.to_numeric(
        metadata["entry_resolution"], errors="coerce"
    )
    invalid_resolution = metadata["entry_resolution"].notna() & ~(
        np.isfinite(metadata["entry_resolution"]) & metadata["entry_resolution"].gt(0)
    )
    if invalid_resolution.any():
        raise ValueError("entry metadata has invalid resolutions")
    candidates = candidates.merge(
        metadata,
        on="entry_pdb_id",
        how="left",
        validate="many_to_one",
    )
    candidates["ligand_contacts_unknown"] = (
        candidates[ligand_contact_columns].isna().any(axis=1)
    )
    candidates["ligand_contact_class"] = np.select(
        [
            candidates["chain_num_contacting_other_ligands"].gt(0).fillna(False),
            candidates["chain_num_contacting_artifacts"].gt(0).fillna(False),
            candidates["chain_num_contacting_ions"].gt(0).fillna(False),
        ],
        [3, 2, 1],
        default=0,
    )
    candidates["num_contacting_ligands"] = candidates[ligand_contact_columns].sum(
        axis=1, min_count=1
    )
    candidates = candidates.sort_values(
        [
            "entry_pdb_id",
            "chain_asym_id",
            "ligand_contacts_unknown",
            "ligand_contact_class",
            "num_contacting_ligands",
            "entry_resolution",
            "biounit_id",
            "chain_instance",
        ],
        ignore_index=True,
        na_position="last",
    ).drop_duplicates(["entry_pdb_id", "chain_asym_id"])
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
            "chain_num_contacting_proteins": "source_num_contacting_proteins",
            "chain_num_contacting_ions": "source_num_contacting_ions",
            "chain_num_contacting_artifacts": "source_num_contacting_artifacts",
            "chain_num_contacting_other_ligands": (
                "source_num_contacting_other_ligands"
            ),
            "entry_resolution": "source_resolution",
        }
    )
    return candidates.loc[:, CANDIDATE_COLUMNS].reset_index(drop=True)


def _parquet_files(path: str | Path) -> list[Path]:
    source = Path(path)
    if source.is_file():
        return [source]
    return sorted(source.rglob("*.parquet")) if source.is_dir() else []


def _sql_string(value: str) -> str:
    return "'" + value.replace("'", "''") + "'"


def _sql_path(path: str | Path) -> str:
    return _sql_string(Path(path).resolve().as_posix())


def _sql_sources(paths: list[Path]) -> str:
    values = ", ".join(_sql_path(path) for path in paths)
    return f"[{values}]"


def _empty_output(path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    pq.write_table(
        pa.Table.from_pylist([], schema=INTERFACE_APO_LINK_SCHEMA),
        temporary,
        compression="zstd",
    )
    temporary.replace(path)
    return path


def write_interface_apo_structure_table(
    mmseqs_alignments: str | Path,
    *,
    foldseek_alignments: str | Path,
    queries: TableInput,
    candidates: TableInput,
    output_path: str | Path,
    config: InterfaceApoSelectionConfig | None = None,
    scratch_dir: str | Path | None = None,
    threads: int = 1,
    memory_limit: str = "7GB",
) -> Path:
    """Write ranked protein-contact-free chain links for interface sides."""
    import duckdb

    if threads < 1:
        raise ValueError("threads must be positive")
    config = config or InterfaceApoSelectionConfig()
    output = Path(output_path)
    query_frame = _read_columns(queries, QUERY_COLUMNS, "interface apo queries")
    candidate_frame = _read_columns(
        candidates, CANDIDATE_COLUMNS, "interface apo candidates"
    )
    if query_frame.empty or candidate_frame.empty:
        return _empty_output(output)
    mmseqs_files = _parquet_files(mmseqs_alignments)
    if not mmseqs_files:
        raise FileNotFoundError(
            f"no MMseqs alignment shards found at {mmseqs_alignments}"
        )
    for path in mmseqs_files:
        missing = sorted(set(ALIGNMENT_COLUMNS).difference(pq.read_schema(path).names))
        if missing:
            raise ValueError(f"MMseqs alignment {path} is missing columns {missing}")
    foldseek_files = _parquet_files(foldseek_alignments)
    for path in foldseek_files:
        required = set(ALIGNMENT_COLUMNS) | {"lddt"}
        missing = sorted(required.difference(pq.read_schema(path).names))
        if missing:
            raise ValueError(f"Foldseek alignment {path} is missing columns {missing}")

    scratch = Path(scratch_dir or output.parent / ".interface-apo-scratch")
    scratch.mkdir(parents=True, exist_ok=True)
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_suffix(output.suffix + ".tmp")
    temporary.unlink(missing_ok=True)
    foldseek_sql = (
        f"""
        SELECT query_entry, target_entry, query_chain_mapped,
               target_chain_mapped, MAX(lddt) AS foldseek_lddt
        FROM read_parquet({_sql_sources(foldseek_files)}, union_by_name=true)
        GROUP BY query_entry, target_entry, query_chain_mapped, target_chain_mapped
        """
        if foldseek_files
        else """
        SELECT NULL::VARCHAR AS query_entry, NULL::VARCHAR AS target_entry,
               NULL::VARCHAR AS query_chain_mapped,
               NULL::VARCHAR AS target_chain_mapped,
               NULL::DOUBLE AS foldseek_lddt
        WHERE FALSE
        """
    )
    connection = duckdb.connect()
    try:
        connection.register("interface_apo_queries", query_frame)
        connection.register("interface_apo_candidates", candidate_frame)
        connection.sql(f"SET threads={threads}")
        connection.sql(f"SET temp_directory={_sql_path(scratch)}")
        connection.sql(f"SET memory_limit={_sql_string(memory_limit)}")
        connection.sql("SET preserve_insertion_order=false")
        connection.sql(
            f"""
            COPY (
                WITH foldseek AS ({foldseek_sql}),
                qualifying_mmseqs AS (
                    SELECT
                        *,
                        ROW_NUMBER() OVER (
                            PARTITION BY query_entry, target_entry,
                                         query_chain_mapped, target_chain_mapped
                            ORDER BY
                                fident DESC,
                                LEAST(qcov, tcov) DESC,
                                qcov DESC,
                                tcov DESC
                        ) AS alignment_rank
                    FROM read_parquet(
                        {_sql_sources(mmseqs_files)}, union_by_name=true
                    )
                    WHERE fident >= {config.min_sequence_identity}
                      AND qcov >= {config.min_query_coverage}
                      AND tcov >= {config.min_target_coverage}
                ),
                mmseqs AS (
                    SELECT
                        query_entry, target_entry, query_chain_mapped,
                        target_chain_mapped,
                        fident AS mmseqs_fident,
                        qcov AS mmseqs_query_coverage,
                        tcov AS mmseqs_target_coverage
                    FROM qualifying_mmseqs
                    WHERE alignment_rank = 1
                ),
                eligible AS (
                    SELECT
                        query.reference_system_id,
                        query.reference_side,
                        query.reference_chain_instance,
                        query.reference_chain_asym_id,
                        candidate.*,
                        mmseqs.mmseqs_fident,
                        mmseqs.mmseqs_query_coverage,
                        mmseqs.mmseqs_target_coverage,
                        foldseek.foldseek_lddt
                    FROM interface_apo_queries AS query
                    INNER JOIN mmseqs
                      ON query.query_entry = mmseqs.query_entry
                     AND query.reference_chain_asym_id
                         = mmseqs.query_chain_mapped
                    INNER JOIN interface_apo_candidates AS candidate
                      ON candidate.source_entry_id = mmseqs.target_entry
                     AND candidate.source_chain_asym_id
                         = mmseqs.target_chain_mapped
                    LEFT JOIN foldseek
                      ON foldseek.query_entry = mmseqs.query_entry
                     AND foldseek.target_entry = mmseqs.target_entry
                     AND foldseek.query_chain_mapped = mmseqs.query_chain_mapped
                     AND foldseek.target_chain_mapped = mmseqs.target_chain_mapped
                    WHERE query.query_entry != candidate.source_entry_id
                      AND candidate.source_num_contacting_proteins = 0
                ),
                ranked AS (
                    SELECT *, ROW_NUMBER() OVER (
                        PARTITION BY reference_system_id, reference_side
                        ORDER BY
                            source_num_contacting_ions IS NULL
                                OR source_num_contacting_artifacts IS NULL
                                OR source_num_contacting_other_ligands IS NULL,
                            CASE
                                WHEN source_num_contacting_other_ligands > 0 THEN 3
                                WHEN source_num_contacting_artifacts > 0 THEN 2
                                WHEN source_num_contacting_ions > 0 THEN 1
                                ELSE 0
                            END,
                            COALESCE(source_num_contacting_ions, 0)
                                + COALESCE(source_num_contacting_artifacts, 0)
                                + COALESCE(
                                    source_num_contacting_other_ligands, 0
                                ),
                            source_resolution IS NULL,
                            source_resolution,
                            mmseqs_fident DESC,
                            LEAST(
                                mmseqs_query_coverage,
                                mmseqs_target_coverage
                            ) DESC,
                            foldseek_lddt DESC NULLS LAST,
                            source_entry_id,
                            source_chain_asym_id
                    ) AS link_rank
                    FROM eligible
                )
                SELECT
                    CAST(reference_system_id AS VARCHAR) AS reference_system_id,
                    CAST(reference_side AS TINYINT) AS reference_side,
                    CAST(reference_chain_instance AS VARCHAR)
                        AS reference_chain_instance,
                    CAST(reference_chain_asym_id AS VARCHAR)
                        AS reference_chain_asym_id,
                    CAST(target_system AS VARCHAR) AS linked_structure_id,
                    CAST(source_entry_id AS VARCHAR) AS source_entry_id,
                    CAST(source_chain_asym_id AS VARCHAR)
                        AS source_chain_asym_id,
                    CAST(source_chain_auth_id AS VARCHAR)
                        AS source_chain_auth_id,
                    CAST(source_biounit_id AS VARCHAR) AS source_biounit_id,
                    CAST(source_chain_instance AS VARCHAR)
                        AS source_chain_instance,
                    CAST(source_num_contacting_proteins AS SMALLINT)
                        AS source_num_contacting_proteins,
                    CAST(source_num_contacting_ions AS SMALLINT)
                        AS source_num_contacting_ions,
                    CAST(source_num_contacting_artifacts AS SMALLINT)
                        AS source_num_contacting_artifacts,
                    CAST(source_num_contacting_other_ligands AS SMALLINT)
                        AS source_num_contacting_other_ligands,
                    CAST(source_resolution AS FLOAT) AS source_resolution,
                    CAST(link_rank AS SMALLINT) AS rank,
                    CAST(mmseqs_fident AS FLOAT) AS mmseqs_fident,
                    CAST(mmseqs_query_coverage AS FLOAT)
                        AS mmseqs_query_coverage,
                    CAST(mmseqs_target_coverage AS FLOAT)
                        AS mmseqs_target_coverage,
                    CAST(foldseek_lddt AS FLOAT) AS foldseek_lddt
                FROM ranked
                WHERE link_rank <= {config.max_per_side}
                ORDER BY reference_system_id, reference_side, rank
            ) TO {_sql_path(temporary)} (
                FORMAT PARQUET, COMPRESSION ZSTD, ROW_GROUP_SIZE 500000
            )
            """
        )
        temporary.replace(output)
    finally:
        connection.close()
        temporary.unlink(missing_ok=True)
    return output


__all__ = [
    "InterfaceApoSelectionConfig",
    "build_interface_apo_candidate_manifest",
    "build_interface_apo_query_manifest",
    "write_interface_apo_structure_table",
]
