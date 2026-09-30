# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
"""Select an unambiguous BindingDB measurement for a receptor."""

from collections.abc import Collection, Mapping, Sequence
from pathlib import Path
from shutil import copyfile
from typing import Any

import numpy as np
import pandas as pd

from plinder.data.annotations.protein_utils import sequences_match_core


def candidate_map(candidates: pd.DataFrame) -> dict[str, list[dict[str, Any]]]:
    """Index strict BindingDB candidates by reported PDB/ligand code."""
    by_key: dict[str, list[dict[str, Any]]] = {}
    for key, target_sequence, endpoint, pchembl, count in candidates[
        ["pdbid_ligid", "target_sequence", "endpoint", "pchembl", "count"]
    ].itertuples(index=False, name=None):
        by_key.setdefault(key, []).append(
            {
                "target_sequence": target_sequence,
                "endpoint": endpoint,
                "pchembl": pchembl,
                "count": count,
            }
        )
    return by_key


def matched_affinity(
    candidates: Sequence[Mapping[str, Any]], receptor_seqres: Mapping[str, str]
) -> tuple[float, str, int] | None:
    """Return pKi/pKd, endpoint, and measurement count for one target match.

    Candidates have already passed the within-target assay/qualifier checks.
    A missing receptor sequence or multiple matching target sequences is
    ambiguous, so neither produces a scalar affinity.
    """
    if not receptor_seqres:
        return None
    matching = [
        candidate
        for candidate in candidates
        if candidate.get("target_sequence")
        and any(
            sequences_match_core(candidate["target_sequence"], sequence)
            for sequence in receptor_seqres.values()
        )
    ]
    if len(matching) != 1:
        return None
    candidate = matching[0]
    return (
        float(candidate["pchembl"]),
        str(candidate["endpoint"]),
        int(candidate["count"]),
    )


def build_ligand_affinity_table(
    annotation_path: Path,
    entry_chains_path: Path,
    candidates_path: Path,
    *,
    pdb_ids: Collection[str] | None = None,
) -> pd.DataFrame:
    """Assign BindingDB candidates to ligands using their receptor chains.

    A BindingDB PDB/CCD cross-reference alone is not sufficient evidence for
    the scalar value. Unmatched or ambiguous ligands retain null values here.
    """
    filters = [("entry_pdb_id", "in", sorted(pdb_ids))] if pdb_ids is not None else None
    annotations = pd.read_parquet(
        annotation_path,
        columns=[
            "ligand_id",
            "system_id",
            "entry_pdb_id",
            "ligand_ccd_code",
            "ligand_protein_chains_asym_id",
        ],
        filters=filters,
    )
    chains = pd.read_parquet(
        entry_chains_path,
        columns=["entry_pdb_id", "chain_asym_id", "chain_sequence"],
        filters=filters,
    )
    sequence_by_chain = {
        (pdb_id, asym_id): sequence
        for pdb_id, asym_id, sequence in chains.itertuples(index=False, name=None)
        if isinstance(sequence, str) and sequence
    }
    candidates_by_key = candidate_map(pd.read_parquet(candidates_path))

    values: list[tuple[float | None, str | None, int | None]] = []
    for _, _, pdb_id, ccd_code, receptor_chains in annotations.itertuples(
        index=False, name=None
    ):
        key = f"{pdb_id}_{ccd_code}".upper()
        receptor_seqres = {}
        if isinstance(receptor_chains, (list, tuple, np.ndarray)):
            for chain in receptor_chains:
                sequence = sequence_by_chain.get(
                    (pdb_id, chain.split(".", maxsplit=1)[-1])
                )
                if sequence:
                    receptor_seqres[chain] = sequence
        values.append(
            matched_affinity(candidates_by_key.get(key, []), receptor_seqres)
            or (None, None, None)
        )
    result = annotations[["ligand_id", "system_id"]].copy()
    affinity_columns = [
        "ligand_binding_affinity",
        "ligand_binding_affinity_endpoint",
        "ligand_binding_affinity_measurement_count",
    ]
    result[affinity_columns] = pd.DataFrame(
        values, index=result.index, columns=affinity_columns
    )
    result["ligand_binding_affinity"] = result["ligand_binding_affinity"].astype(float)
    result["ligand_binding_affinity_endpoint"] = result[
        "ligand_binding_affinity_endpoint"
    ].astype("string")
    result["ligand_binding_affinity_measurement_count"] = result[
        "ligand_binding_affinity_measurement_count"
    ].astype("Int64")
    result["system_has_binding_affinity"] = (
        result.groupby("system_id")["ligand_binding_affinity"].transform("count").gt(0)
    )
    return result.drop(columns="system_id")


def publish_affinity_tables(
    data_dir: Path,
    *,
    affinity_dir: Path | None = None,
    base_dir: Path | None = None,
    changed_pdb_ids: Collection[str] = (),
    removed_pdb_ids: Collection[str] = (),
) -> None:
    """Write corrected ligand values and source measurements beside the index."""
    affinity_dir = affinity_dir or data_dir / "dbs" / "affinity"
    index_dir = data_dir / "index"
    index_dir.mkdir(parents=True, exist_ok=True)
    ligand_path = index_dir / "ligand_affinity.parquet"
    records_path = index_dir / "bindingdb_measurements.parquet"
    ligand_tmp = ligand_path.with_suffix(".parquet.tmp")
    records_tmp = records_path.with_suffix(".parquet.tmp")
    try:
        if base_dir is None:
            values = build_ligand_affinity_table(
                index_dir / "annotation_table.parquet",
                index_dir / "entry_chains.parquet",
                affinity_dir / "candidates.parquet",
            )
        else:
            old = pd.read_parquet(base_dir / "index/ligand_affinity.parquet")
            if removed_pdb_ids:
                removed_ligands = pd.read_parquet(
                    base_dir / "index/annotation_table.parquet",
                    columns=["ligand_id"],
                    filters=[("entry_pdb_id", "in", sorted(removed_pdb_ids))],
                )["ligand_id"]
                old = old.loc[~old["ligand_id"].isin(removed_ligands)]
            new = (
                build_ligand_affinity_table(
                    index_dir / "annotation_table.parquet",
                    index_dir / "entry_chains.parquet",
                    affinity_dir / "candidates.parquet",
                    pdb_ids=changed_pdb_ids,
                )
                if changed_pdb_ids
                else old.iloc[:0]
            )
            ligand_ids = pd.read_parquet(
                index_dir / "annotation_table.parquet", columns=["ligand_id"]
            )["ligand_id"]
            values = (
                pd.concat([old, new], ignore_index=True)
                .set_index("ligand_id", verify_integrity=True)
                .loc[ligand_ids]
                .reset_index()
            )
        values.to_parquet(ligand_tmp, index=False)
        copyfile(affinity_dir / "measurements.parquet", records_tmp)
        ligand_tmp.replace(ligand_path)
        records_tmp.replace(records_path)
    finally:
        ligand_tmp.unlink(missing_ok=True)
        records_tmp.unlink(missing_ok=True)
