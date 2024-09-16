# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from biotite import structure as struc
from biotite.structure import AtomArray

from plinder.core.structure.atoms import (
    atom_array_from_cif_file,
)
from plinder.core.structure.contacts import (
    get_atom_neighbors,
    label_structure_gaps,
    pairwise_chain_contacts,
)
from plinder.core.utils.log import setup_logger

log = setup_logger(__name__)


def get_contacts_gaps_overlap(
    atoms: AtomArray,
    gaps: pd.DataFrame | None,
    contacts: pd.DataFrame,
    pdb_id: str,
    smaller_radius: float = 4.0,
    larger_radius: float = 8.0,
) -> pd.DataFrame | None:
    smaller_radius = 4.0
    larger_radius = 8.0
    annotations = []

    for i, r in contacts.iterrows():
        ch1 = r.chain1
        ch2 = r.chain2
        resid1 = r.resid1
        resid2 = r.resid2
        conts = r.contacts

        annotation = {
            "pdb_id": pdb_id,
            "chain1": ch1,
            "resid1": resid1,
            "chain2": ch2,
            "resid2": resid2,
        }
        if len(conts) == 0:
            continue

        ch1_resi = set([p[0] for p in conts])
        ch2_resi = set([p[1] for p in conts])

        interface_mask = (
            (atoms.chain_id == ch1) & np.isin(atoms.res_id, list(ch1_resi))
        ) | ((atoms.chain_id == ch2) & np.isin(atoms.res_id, list(ch2_resi)))
        interface = atoms[interface_mask].copy()
        if gaps is None:
            continue
        ch1_gaps = gaps.query(f'chain == "{ch1}"')
        ch2_gaps = gaps.query(f'chain == "{ch2}"')
        ch1_gaps = set(ch1_gaps.gap_start).union(set(ch1_gaps.gap_end))
        ch2_gaps = set(ch2_gaps.gap_start).union(set(ch2_gaps.gap_end))
        gap_mask = ((atoms.chain_id == ch1) & np.isin(atoms.res_id, list(ch1_gaps))) | (
            (atoms.chain_id == ch2) & np.isin(atoms.res_id, list(ch2_gaps))
        )
        gap_atoms = atoms[gap_mask].copy()

        if not gap_atoms.shape[0]:
            continue

        for radius in [smaller_radius, larger_radius]:
            rad_label = int(round(radius, 0))
            gap_neigh = get_atom_neighbors(interface, gap_atoms, radius=radius)
            annot = pd.DataFrame(
                [
                    {"resi": at.res_id, "resn": at.res_name, "chain": at.chain_id}
                    for at in gap_neigh
                ]
            )
            missing_interface_res = 0
            if annot.empty:
                annotation[f"interface_atom_gaps_{rad_label}A"] = 0
                annotation[f"missing_interface_residues_{rad_label}A"] = 0
            else:
                for ch, ch_gap in [(ch1, ch1_gaps), (ch2, ch2_gaps)]:
                    missing_interface_res += len(
                        set(annot.query(f'chain == "{ch}"').resi).intersection(ch_gap)
                    )
                annotation[f"interface_atom_gaps_{rad_label}A"] = annot.shape[0]
                annotation[
                    f"missing_interface_residues_{rad_label}A"
                ] = missing_interface_res

        annotations.append(annotation)

    return annotations


# TODO: review this function
# it does not use the ligand chain definitions!
def annotate_interface_gaps(
    cif_file: Path,
    protein_chains: list[str] | None = None,
    ligand_chains: list[str] | None = None,
    smaller_radius: float = 4.0,
    larger_radius: float = 8.0,
    use_author_fields: bool = False,
) -> dict[str, Any]:
    pdb_id = cif_file.stem.split("_")[0].replace("pdb_0000", "")

    if ".cif" in cif_file.name:
        atoms = atom_array_from_cif_file(cif_file, use_author_fields=use_author_fields)
    else:
        raise ValueError(f"unsupported file extension: {cif_file}")
    assert atoms is not None

    # Complex atom array
    lig_filter = atoms.hetero
    prot_filter = struc.filter_amino_acids(atoms)

    if ligand_chains is not None:
        # Filter atoms of interest
        lig_filter = atoms.hetero & np.isin(atoms.chain_id, np.array(ligand_chains))
    if protein_chains is not None:
        prot_filter = struc.filter_amino_acids(atoms) & np.isin(
            atoms.chain_id,
            np.array(protein_chains),
        )
    prot_arr = atoms[prot_filter].copy()
    complex_arr = atoms[prot_filter | lig_filter].copy()
    ppi_contacts, pli_contacts = pairwise_chain_contacts(complex_arr)
    gaps = label_structure_gaps(prot_arr)
    pli_missing_df = pd.DataFrame(
        get_contacts_gaps_overlap(
            atoms,
            gaps,
            pli_contacts,
            pdb_id,
            smaller_radius,
            larger_radius,
        )
    )
    ppi_missing_df = pd.DataFrame(
        get_contacts_gaps_overlap(
            atoms,
            gaps,
            ppi_contacts,
            pdb_id,
            smaller_radius,
            larger_radius,
        )
    )
    ppi_missing_dict: dict[str, dict[tuple[str, str], dict[str, int]]] = {}
    pli_missing_dict: dict[str, dict[tuple[str, str], dict[str, int]]] = {}
    if not ppi_missing_df.empty:
        ppi_missing_df["chain1_chain2"] = ppi_missing_df[["chain1", "chain2"]].apply(
            tuple, axis=1
        )
        ppi_missing_df.set_index("chain1_chain2", inplace=True)
        if ppi_missing_df.index.duplicated().any():
            log.warn(f"found duplicate chain1_chain2 in ppi_missing_df for {cif_file}")
            ppi_missing_df = ppi_missing_df[
                ~ppi_missing_df.index.duplicated(keep="first")
            ]
        ppi_missing_dict = ppi_missing_df[
            [
                "interface_atom_gaps_4A",
                "missing_interface_residues_4A",
                "interface_atom_gaps_8A",
                "missing_interface_residues_8A",
            ]
        ].to_dict(orient="index")
    if not pli_missing_df.empty:
        pli_missing_df["chain1_chain2"] = pli_missing_df[["chain1", "chain2"]].apply(
            tuple, axis=1
        )

        pli_missing_df.set_index("chain1_chain2", inplace=True)
        if pli_missing_df.index.duplicated().any():
            log.warn(f"found duplicate chain1_chain2 in pli_missing_df for {cif_file}")
            pli_missing_df = pli_missing_df[
                ~pli_missing_df.index.duplicated(keep="first")
            ]
        pli_missing_dict = pli_missing_df[
            [
                "interface_atom_gaps_4A",
                "missing_interface_residues_4A",
                "interface_atom_gaps_8A",
                "missing_interface_residues_8A",
            ]
        ].to_dict(orient="index")

    return {
        "ppi_interface_gap_annotation": ppi_missing_dict,
        "ligand_interface_gap_annotation": pli_missing_dict,
    }
