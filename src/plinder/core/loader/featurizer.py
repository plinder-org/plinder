import random
from typing import Any

import torch

from plinder.core.index.system import PlinderSystem
from plinder.core.loader.utils import pad_and_stack
from plinder.core.structure.atoms import (
    _one_hot_encode_stack,
    _stack_ligand_feat,
)
from plinder.core.structure.diffdock_utils import lig_atom_featurizer
from plinder.core.utils import constants as pc


def system_featurizer(
    system: PlinderSystem,
    pad_value: int = -100,
    featurize_apo: bool = True,
    seed: int = 42,
) -> dict[str, Any]:
    # Set seed
    random.seed(seed)
    # Load holo and alternate (apo and pred) structures
    holo_structure = system.holo_structure
    apo_structures = system.alternate_structures

    # This must be used to order the chain features
    protein_chain_order = holo_structure.protein_chain_ordered
    ligand_chain_order = holo_structure.ligand_chain_ordered

    # Extract holo chain id
    assert holo_structure.protein_sequence is not None
    holo_chain = list(holo_structure.protein_sequence.keys())[0]

    # Get input sequence dictionary
    input_sequences = holo_structure.protein_sequence
    # Get residue-level features and sort in input structure order
    input_sequence_residue_feat_stack = _one_hot_encode_stack(
        [input_sequences[ch] for ch in protein_chain_order], pc.AA_TO_INDEX, "UNK"
    )
    # Get atom-level features and sort in input structure order
    input_sequence_full_atom_dict = holo_structure.input_sequence_full_atom_dict
    input_sequence_full_atom_feat_stack = _one_hot_encode_stack(
        [input_sequence_full_atom_dict[ch] for ch in protein_chain_order],
        pc.ELE2NUM,
        "other",
    )

    sequence_features = {
        "input_sequence_residue_feat_stack": input_sequence_residue_feat_stack,
        "input_sequence_full_atom_feat_stack": input_sequence_full_atom_feat_stack,
    }

    # Get holo structure to input sequence atom
    holo_sequence_atom_mask_stacked = holo_structure.sequence_atom_mask
    # Get holo structure to input sequence atom
    holo_input_sequence_residue_mask_stacked = (
        holo_structure.input_sequence_residue_mask_stacked
    )
    # Get holo structure coordinates
    holo_protein_coordinates_stacked = holo_structure.protein_coords
    # Get holo calpha coordinates
    holo_protein_calpha_coordinates_stacked = holo_structure.protein_calpha_coords

    # Bundle holo features into a dictionary
    holo_features = {
        "holo_sequence_atom_mask_stacked": holo_sequence_atom_mask_stacked,
        "holo_input_sequence_residue_mask_stacked": holo_input_sequence_residue_mask_stacked,
        "holo_protein_coordinates_stacked": holo_protein_coordinates_stacked,
        "holo_protein_calpha_coordinates_stacked": holo_protein_calpha_coordinates_stacked,
    }
    if featurize_apo:
        selected_apo = apo_structures[random.choice(list(apo_structures.keys()))]
        # Set apo chain to match holo
        selected_apo.set_chain(holo_chain)
        apo_sequence_atom_mask_stacked = selected_apo.sequence_atom_mask
        apo_input_sequence_residue_mask_stacked = (
            selected_apo.input_sequence_residue_mask_stacked
        )
        apo_protein_coordinates_array_stacked = selected_apo.protein_coords
        apo_protein_calpha_coordinates_array_stacked = (
            selected_apo.protein_calpha_coords
        )
        # Apo to holo cropping mask
        apo_sequence_to_holo_structure_mask_stacked = (
            holo_structure.protein_structure_residue_mask(selected_apo)
        )

        apo_features = {
            "apo_sequence_atom_masks_stacked": apo_sequence_atom_mask_stacked,
            "apo_input_sequence_residue_masks_stacked": apo_input_sequence_residue_mask_stacked,
            "apo_protein_coordinates_arrays_stacked": apo_protein_coordinates_array_stacked,
            "apo_protein_calpha_coordinates_arrays_stacked": apo_protein_calpha_coordinates_array_stacked,
            "apo_sequence_to_holo_structure_mask_stacked": apo_sequence_to_holo_structure_mask_stacked,
        }
    else:
        apo_features = {
            "apo_sequence_atom_masks_stacked": [],
            "apo_input_sequence_residue_masks_stacked": [],
            "apo_protein_coordinates_arrays_stacked": [],
            "apo_protein_calpha_coordinates_arrays_stacked": [],
            "apo_sequence_to_holo_structure_mask_stacked": [],
        }

    # Ligand features
    input_ligand_templates = (
        holo_structure.input_ligand_templates
    )  # 2D templates from SMILES

    input_ligand_conformers_coords = (
        holo_structure.input_ligand_conformers_coords
    )  # 3D (random) conformer
    resolved_ligand_mols_coords = holo_structure.resolved_ligand_mols_coords

    input_ligand_feat = {
        ch: lig_atom_featurizer(ligand_mol)
        for ch, ligand_mol in input_ligand_templates.items()
    }
    # Stack in ligand_chain_order order
    input_ligand_feat_stack = [
        feats for feats in _stack_ligand_feat(input_ligand_feat, ligand_chain_order)
    ]

    input_conformer_ligand_coords_stack = [
        coord
        for coord in _stack_ligand_feat(
            input_ligand_conformers_coords, ligand_chain_order
        )
    ]

    # Get resolved ligand mols coordinate
    resolved_ligand_mols_coords_stack = [
        coord
        for coord in _stack_ligand_feat(resolved_ligand_mols_coords, ligand_chain_order)
    ]
    ligand_features = {
        "input_conformer_ligand_feat_stack": input_ligand_feat_stack,
        "input_conformer_ligand_coords_stack": input_conformer_ligand_coords_stack,
        "resolved_ligand_mols_coords_stack": resolved_ligand_mols_coords_stack,
    }

    features = {
        "sequence_features": sequence_features,
        "holo_features": holo_features,
        "apo_features": apo_features,
        "ligand_features": ligand_features,
    }

    # Pad tensors to make chains have equal length
    padded_features = {
        feat_group_name: {
            feat_name: (
                pad_and_stack(
                    [torch.tensor(feat_per_chain) for feat_per_chain in feat],
                    dim=0,
                    value=pad_value,
                )
                if len(feat) != 0
                else []
            )
            for feat_name, feat in feat_dict.items()
        }
        for feat_group_name, feat_dict in features.items()
    }

    # Set features as new properties
    return padded_features
