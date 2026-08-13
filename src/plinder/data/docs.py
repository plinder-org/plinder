# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
from __future__ import annotations

from functools import cache
from pathlib import Path

import pandas as pd
import pyarrow.parquet as pq

from plinder.core.release import RELEASE_TABLES, PlinderRelease
from plinder.data import column_descriptions

TSV_DIR = Path(column_descriptions.__file__).parent
TABLE_TSV_DIR = TSV_DIR / "tables"

DERIVED_COLUMN_DESCRIPTIONS = {
    "chain_is_ligand_like": (
        "Whether this chain is ligand-like rather than an eligible receptor chain"
    ),
    "chain_type": "Raw CIF receptor polymer type",
    "chain_receptor_type": ("Receptor polymer category: protein, dna, rna, or dna+rna"),
    "chain_is_holo": (
        "Whether the chain participates in a holo receptor or is excluded from "
        "the apo partition"
    ),
    "chain_uniprot_ids": "UniProt accessions mapped to the chain",
    "chain_sequence": (
        "Full polymer sequence in one-letter code, keyed by the source "
        "asymmetric-chain ID"
    ),
    "biounit_id": "Biological assembly identifier",
    "chain_instance": (
        "Assembly chain instance encoded as <operation>.<label_asym_id>"
    ),
    "chain_role": "Assembly role: receptor, ligand, or water",
    "chain_num_contacting_ions": (
        "Number of ion chains contacting this assembly-chain instance"
    ),
    "chain_num_contacting_artifacts": (
        "Number of crystallization-artifact chains contacting this "
        "assembly-chain instance"
    ),
    "chain_num_contacting_other_ligands": (
        "Number of other ligand chains contacting this assembly-chain instance"
    ),
    "source_mmcif_major_revision": (
        "Major revision number of the source PDB mmCIF used during ingest"
    ),
    "source_mmcif_minor_revision": (
        "Minor revision number of the source PDB mmCIF used during ingest"
    ),
    "reference_system_id": "Holo PLINDER system linked to this apo chain",
    "linked_structure_id": "Stable identifier for the linked apo chain",
    "source_entry_id": "PDB entry containing the linked apo chain",
    "source_chain_asym_id": "Source asymmetric-chain ID of the linked apo chain",
    "source_chain_auth_id": "Author chain ID of the linked apo chain",
    "source_biounit_id": "Source biological assembly containing the apo chain",
    "source_chain_instance": "Exact biological-assembly chain instance that was scored",
    "source_num_contacting_ions": "Number of ion chains contacting the apo chain",
    "source_num_contacting_artifacts": (
        "Number of crystallization-artifact chains contacting the apo chain"
    ),
    "source_num_contacting_other_ligands": (
        "Number of other ligand chains contacting the apo chain"
    ),
    "source_resolution": "Experimental resolution of the apo source entry",
    "rank": "Apo candidate rank within the reference holo system",
    "num_ligand_pockets": "Number of holo ligand pockets matched by this apo chain",
    "min_pocket_fident": "Minimum pocket sequence identity across matched ligands",
    "mean_pocket_fident": "Mean pocket sequence identity across matched ligands",
    "min_protein_fident_weighted_sum": (
        "Minimum chain-length-weighted protein sequence identity across matched ligands"
    ),
    "min_protein_fident_qcov_weighted_sum": (
        "Minimum chain-length-weighted protein identity times query coverage "
        "across matched ligands"
    ),
    "min_protein_lddt_weighted_sum": (
        "Minimum chain-length-weighted protein LDDT across matched ligands"
    ),
    "selected_residue_numbers": (
        "Source-mmCIF residue numbers selected by any ligand pocket or protein "
        "interface, ordered by selected_residue_indices"
    ),
    "selected_residue_indices": (
        "Zero-based resolved-chain indices corresponding positionally to "
        "selected_residue_numbers"
    ),
    "interface_chain_1": "First canonical assembly-chain instance",
    "interface_chain_2": "Second canonical assembly-chain instance",
    "interface_chain_1_residue_numbers": (
        "Source-mmCIF residue numbers at the interface on chain 1"
    ),
    "interface_chain_1_residue_indices": (
        "Zero-based resolved-chain indices corresponding positionally to chain "
        "1 residue numbers"
    ),
    "interface_chain_2_residue_numbers": (
        "Source-mmCIF residue numbers at the interface on chain 2"
    ),
    "interface_chain_2_residue_indices": (
        "Zero-based resolved-chain indices corresponding positionally to chain "
        "2 residue numbers"
    ),
    "interface_num_contact_residue_pairs": (
        "Number of residue pairs in contact across the protein interface"
    ),
    "prodigy_is_annotated": (
        "Whether PRODIGY-cryst features and a biological/crystal label are available"
    ),
    "prodigy_label": "PRODIGY-cryst BIO or XTAL interface label",
    "prodigy_probability_bio": (
        "PRODIGY-cryst probability that the interface is biological"
    ),
    "prodigy_link_density": "PRODIGY-cryst interface link-density feature",
    "prodigy_intermolecular_contacts": (
        "PRODIGY-cryst count of intermolecular atomic contacts"
    ),
    "prodigy_charged_charged_contacts": (
        "PRODIGY-cryst count of charged-charged contacts"
    ),
    "prodigy_charged_polar_contacts": ("PRODIGY-cryst count of charged-polar contacts"),
    "prodigy_charged_apolar_contacts": (
        "PRODIGY-cryst count of charged-apolar contacts"
    ),
    "prodigy_polar_polar_contacts": ("PRODIGY-cryst count of polar-polar contacts"),
    "prodigy_apolar_polar_contacts": ("PRODIGY-cryst count of apolar-polar contacts"),
    "prodigy_apolar_apolar_contacts": ("PRODIGY-cryst count of apolar-apolar contacts"),
    "representative_ligand_id": (
        "Ligand ID selected as the representative for a ligand-pocket group"
    ),
    "representative_system_id": (
        "System ID selected as the representative for a ligand pocket or "
        "protein interface"
    ),
    "receptor_chain_asym_ids": (
        "Asymmetric-chain IDs in the representative receptor-chain set"
    ),
    "receptor_set_id": (
        "Stable identifier for the representative set of receptor chains"
    ),
    "pocket_residues": (
        "Representative pocket residues encoded as asym ID, residue number, "
        "and zero-based resolved-chain index"
    ),
    "interactions": (
        "Representative protein-ligand interactions encoded with their receptor "
        "residue mapping and interaction type"
    ),
    "half_interface_id": "Stable identifier for one representative interface side",
    "instance_chain_id": "Assembly-chain instance forming this half-interface",
    "residue_numbers": "Source-mmCIF residue numbers in this half-interface",
    "residue_indices": (
        "Zero-based resolved-chain indices corresponding positionally to "
        "residue_numbers"
    ),
    "half_interface_1_id": "First half-interface of this representative",
    "half_interface_2_id": "Second half-interface of this representative",
    "side_1_half_interface_id": (
        "Representative half-interface assigned to the first query side"
    ),
    "side_2_half_interface_id": (
        "Representative half-interface assigned to the second query side"
    ),
    "uniqueness": (
        "Identifier differentiating systems that are simple crystal symmetries "
        "within a biological assembly"
    ),
    "biounit_num_ligands": "Number of ligands in the biological assembly",
    "biounit_num_unique_ccd_codes": (
        "Number of distinct ligand CCD codes in the biological assembly"
    ),
    "biounit_num_proper_ligands": (
        "Number of proper ligands in the biological assembly"
    ),
    "system_protein_chains_total_length": (
        "Total length of all receptor polymer chains in the system"
    ),
    "system_unique_ccd_codes": "Distinct ligand CCD codes in the system",
    "system_proper_unique_ccd_codes": (
        "Distinct CCD codes of proper ligands in the system"
    ),
    "ligand_is_3d_score_able": (
        "Whether the canonical ligand SDF supports finite shape, color, and "
        "SuCOS self-scoring"
    ),
    "ligand_smiles_id": (
        "Integer node ID assigned to this exact canonical SMILES for ligand "
        "similarity scoring"
    ),
    "ligand_max_cofactor_similarity": (
        "Maximum ECFP4/1024 Tanimoto similarity, as a percentage, to any CCD "
        "structure in the cofactor list"
    ),
    "ligand_most_similar_cofactor": (
        "CCD code of the cofactor with maximum ECFP4/1024 Tanimoto similarity"
    ),
    "ligand_is_cofactor_like": (
        "Whether maximum similarity to a listed CCD cofactor is at least 90 percent"
    ),
    "ligand_tanimoto_ecfp4_1024_90_cluster": (
        "Cluster ID from the 90-percent ECFP4/1024 Tanimoto set cover; every "
        "member has a direct threshold-qualified edge to its representative"
    ),
    "ligand_tanimoto_ecfp4_1024_90_cluster_num_pdb_ids": (
        "Number of distinct PDB entries represented in the ligand's 90-percent "
        "Tanimoto set-cover cluster"
    ),
}

DERIVED_COLUMN_DESCRIPTIONS.update(
    {
        f"system_ligand_has_{name}": description
        for name, description in {
            "lipinski": "Whether the system has a Lipinski ligand",
            "cofactor": "Whether the system has a cofactor ligand",
            "fragment": "Whether the system has a fragment ligand",
            "monosaccharide": (
                "Whether the system has a ligand containing one saccharide unit"
            ),
            "oligosaccharide": (
                "Whether the system has a ligand containing multiple saccharide units"
            ),
            "mononucleotide": (
                "Whether the system has a ligand containing one nucleotide unit"
            ),
            "oligonucleotide": (
                "Whether the system has a ligand containing multiple nucleotide units"
            ),
            "monopeptide": (
                "Whether the system has a ligand containing one peptide unit"
            ),
            "oligopeptide": (
                "Whether the system has a ligand containing multiple peptide units"
            ),
            "artifact": "Whether the system has an artifact ligand",
            "other": "Whether the system has a ligand classified as other",
            "covalent": "Whether the system has a covalent ligand",
            "invalid": "Whether the system has an invalid ligand",
            "ion": "Whether the system has an ion",
        }.items()
    }
)


def get_cluster_column_descriptions(
    plindex: pd.DataFrame,
) -> list[tuple[str, str | None, str | None]]:
    rows: list[tuple[str, str | None, str | None]] = []
    set_cover_columns = [
        c
        for c in plindex.columns
        if c.endswith("__set_cover") and not c.endswith("__directed_set_cover")
    ]
    for column in set_cover_columns:
        parts = column.split("__")
        metric, threshold = parts[:2]
        ligand_level = parts[2] == "ligand"
        level = "ligand-level " if ligand_level else ""
        rows.append(
            (
                column,
                "str",
                f"Cluster ID for {level}set cover built from reciprocal-minimum "
                f"{metric} with {threshold} threshold; each member has a direct "
                "threshold-qualified edge to its representative",
            )
        )
    directed_cover_columns = [
        c for c in plindex.columns if c.endswith("directed_set_cover")
    ]
    for column in directed_cover_columns:
        parts = column.split("__")
        metric, threshold = parts[:2]
        ligand_level = parts[2] == "ligand"
        half_interface = parts[-1].startswith("chain_")
        level = (
            "ligand-level "
            if ligand_level
            else f"{parts[-1].removesuffix('_directed_set_cover').replace('_', ' ')} "
            if half_interface
            else ""
        )
        rows.append(
            (
                column,
                "str",
                f"Cluster ID for {level}directed set cover built from "
                f"{metric} with {threshold} threshold; each member's "
                "query-to-centroid score meets the threshold",
            )
        )
    centroid_columns = [
        c
        for c in plindex.columns
        if c.endswith(
            (
                "__set_cover__is_centroid",
                "__directed_set_cover__is_centroid",
            )
        )
    ]
    for column in centroid_columns:
        parts = column.split("__")
        metric, threshold = parts[:2]
        ligand_level = parts[2] == "ligand"
        level = "ligand-level " if ligand_level else ""
        cover_kind = (
            "directed set-cover"
            if "__directed_set_cover__" in column
            else "set-cover"
        )
        rows.append(
            (
                column,
                "bool | None",
                f"Whether this row is the published centroid for its {level}"
                f"{cover_kind} cluster built from {metric} with "
                f"{threshold} threshold; missing means the row is outside the "
                "clustering universe",
            )
        )
    directed_coverage_columns = [
        c
        for c in plindex.columns
        if c.endswith(
            (
                "__directed_set_cover__coverage_count",
                "__directed_set_cover__coverage_fraction",
            )
        )
    ]
    for column in directed_coverage_columns:
        parts = column.split("__")
        metric, threshold = parts[:2]
        is_count = column.endswith("__coverage_count")
        quantity = (
            "Number of directed-cover query nodes this ligand could initially cover"
            if is_count
            else "Fraction of its directed weak component this ligand could "
            "initially cover"
        )
        rows.append(
            (
                column,
                "int | None" if is_count else "float | None",
                f"{quantity} at {metric} {threshold} threshold, including "
                "itself; missing means the row is outside the clustering "
                "universe",
            )
        )
    column_order = {column: index for index, column in enumerate(plindex.columns)}
    return sorted(rows, key=lambda row: column_order[row[0]])


@cache
def _model_description_lookup() -> dict[str, str]:
    """Generate descriptions from the annotation models that emit the columns."""
    from plinder.data.annotations.aggregate_annotations import Entry, System
    from plinder.data.annotations.get_ligand_validation import (
        EntryValidation,
        ResidueListValidation,
    )
    from plinder.data.annotations.ligand_utils import Ligand
    from plinder.data.annotations.protein_utils import Chain

    descriptions: dict[str, str] = {}
    model_prefixes = [
        (Entry, "entry"),
        (EntryValidation, "entry_validation"),
        (System, "system"),
        (Ligand, "ligand"),
        (Chain, "chain"),
    ]
    model_prefixes.extend(
        (Chain, prefix)
        for prefix in [
            "system_protein_chains",
            "system_ligand_chains",
            "ligand_interacting_ligand_chains",
            "ligand_neighboring_ligand_chains",
            "ligand_protein_chains",
        ]
    )
    model_prefixes.extend(
        (ResidueListValidation, f"{prefix}_validation")
        for prefix in [
            "system_pocket",
            "system_ligand",
            "ligand_interacting_ligand_chains",
            "ligand_neighboring_ligand_chains",
            "ligand_protein_chains",
            "system_ligand_chains",
            "system_protein_chains",
        ]
    )
    for model, prefix in model_prefixes:
        for name, _, description in model.document_properties(prefix):
            descriptions[name] = description
    return descriptions


def _base_description_lookup() -> dict[str, str]:
    """Combine model-generated and explicitly derived column descriptions."""
    descriptions: dict[str, str] = {}
    descriptions.update(_model_description_lookup())
    descriptions.update(DERIVED_COLUMN_DESCRIPTIONS)
    return descriptions


def get_table_column_descriptions(*, table_name: str, schema) -> pd.DataFrame:
    """Build ordered descriptions for exactly the columns in one release table.

    The Arrow schema supplies the published names, order, and data types.  A
    missing prose description is an error so documentation cannot silently lag
    behind a release table.
    """
    if table_name not in RELEASE_TABLES:
        choices = ", ".join(sorted(RELEASE_TABLES))
        raise KeyError(f"unknown release table {table_name!r}; choose from: {choices}")

    fields = list(schema)
    names = [field.name for field in fields]
    descriptions = _base_description_lookup()
    cluster_rows = get_cluster_column_descriptions(pd.DataFrame(columns=names))
    descriptions.update({name: description for name, _, description in cluster_rows})
    missing = [name for name in names if name not in descriptions]
    if missing:
        raise ValueError(
            f"release table {table_name!r} has undocumented columns: {missing}"
        )
    return pd.DataFrame(
        {
            "Name": names,
            "Type": [str(field.type) for field in fields],
            "Description": [descriptions[name] for name in names],
        }
    )


def get_column_descriptions(
    table_name: str,
    *,
    description_dir: Path = TABLE_TSV_DIR,
) -> pd.DataFrame:
    """Read the checked-in descriptions for one release parquet table."""
    if table_name not in RELEASE_TABLES:
        choices = ", ".join(sorted(RELEASE_TABLES))
        raise KeyError(f"unknown release table {table_name!r}; choose from: {choices}")
    path = description_dir / f"{table_name}.tsv"
    if not path.is_file():
        raise FileNotFoundError(f"missing release column descriptions: {path}")
    frame = pd.read_csv(path, sep="\t")
    if frame["Name"].duplicated().any():
        duplicates = frame.loc[frame["Name"].duplicated(), "Name"].tolist()
        raise ValueError(f"duplicate descriptions in {path}: {duplicates}")
    return frame


def write_column_descriptions(
    *,
    data_dir: Path,
    output_dir: Path = TABLE_TSV_DIR,
) -> None:
    """Write one complete description TSV per published parquet table."""
    release = PlinderRelease(data_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    expected: set[Path] = set()
    for table_name, table in RELEASE_TABLES.items():
        table_path = release.fetch(str(table["artifact"]))
        descriptions = get_table_column_descriptions(
            table_name=table_name,
            schema=pq.read_schema(table_path),
        )
        output_path = output_dir / f"{table_name}.tsv"
        descriptions.to_csv(output_path, sep="\t", index=False)
        expected.add(output_path)
    for stale in output_dir.glob("*.tsv"):
        if stale not in expected:
            stale.unlink()
