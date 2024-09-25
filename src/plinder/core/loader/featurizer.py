from typing import Any

import torch

from plinder.core.loader.utils import pad_and_stack
from plinder.core.structure.atoms import (
    _one_hot_encode_stack,
    _stack_atom_array_features,
    _stack_ligand_feat,
)
from plinder.core.structure.diffdock_utils import lig_atom_featurizer
from plinder.core.structure.structure import Structure
from plinder.core.utils import constants as pc


def structure_featurizer(structure: Structure, pad_value: int = -100) -> dict[str, Any]:
    # This must be used to order the chain features
    protein_chain_order = structure.protein_chain_ordered
    ligand_chain_order = structure.ligand_chain_ordered
    protein_atom_array = structure.protein_atom_array
    sequence_atom_mask_stacked = structure.sequence_atom_mask
    input_sequence_residue_mask_stacked = structure.input_sequence_residue_mask_stacked
    protein_coordinates_stacked = structure.protein_coords
    protein_calpha_coordinates_stacked = structure.protein_calpha_coords
    input_ligand_templates = (
        structure.input_ligand_templates
    )  # 2D templates from SMILES
    # input_ligand_conformers = structure.input_ligand_conformers  # not needed as features are assigned from 2D
    input_ligand_conformers_coords = (
        structure.input_ligand_conformers_coords
    )  # 3D (random) conformer
    resolved_ligand_mols_coords = structure.resolved_ligand_mols_coords

    # Sequence atom-level features
    input_sequence_full_atom_feat = structure.input_sequence_full_atom_feat

    # Featurize and sort in input structure order
    input_sequence_full_atom_feat_stack = _one_hot_encode_stack(
        [input_sequence_full_atom_feat[ch] for ch in protein_chain_order],
        pc.ELE2NUM,
        "other",
    )

    # Get residue type feature
    protein_structure_residue_type_arr = _stack_atom_array_features(
        protein_atom_array, "res_name", protein_chain_order
    )
    protein_structure_residue_type_stack = [
        feat
        for feat in _one_hot_encode_stack(
            protein_structure_residue_type_arr, pc.AA_TO_INDEX, "UNK"
        )
    ]
    # TODO: Fix issues with ligands conformer generation
    # Featurize and stack ligand chains
    # VO: try passing the 2D - does not need a conformer!
    input_conformer_ligand_feat = {
        ch: lig_atom_featurizer(ligand_mol)
        for ch, ligand_mol in input_ligand_templates.items()
    }
    # Stack in ligand_chain_order order
    input_conformer_ligand_feat_stack = [
        feats
        for feats in _stack_ligand_feat(input_conformer_ligand_feat, ligand_chain_order)
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
    features = {
        "sequence_atom_mask_feature": sequence_atom_mask_stacked,
        "input_sequence_residue_mask_feature": input_sequence_residue_mask_stacked,
        "protein_coordinates": protein_coordinates_stacked,
        "protein_calpha_coordinates": protein_calpha_coordinates_stacked,
        "input_sequence_full_atom_feature": input_sequence_full_atom_feat_stack,
        "protein_structure_residue_feature": protein_structure_residue_type_stack,
        "input_conformer_ligand_feature": input_conformer_ligand_feat_stack,
        "input_conformer_ligand_coordinates": input_conformer_ligand_coords_stack,
        "resolved_ligand_mols_feature": resolved_ligand_mols_coords_stack,
    }

    # Pad tensors to make chains have equal length
    padded_features = {
        feat_name: pad_and_stack(
            [torch.tensor(feat_per_chain) for feat_per_chain in feat],
            dim=0,
            value=pad_value,
        )
        for feat_name, feat in features.items()
    }

    # Set features as new properties
    return padded_features
