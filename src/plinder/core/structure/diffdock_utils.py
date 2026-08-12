# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
# mypy: disable-error-code="no-untyped-def, no-untyped-call, attr-defined, assignment, arg-type, var-annotated"
# ruff: noqa

from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")

"""
    Ligand atom featurization taken/adapted from https://github.com/gcorso/DiffDock.

    Only ``lig_atom_featurizer`` is wired into plinder (via
    ``plinder.core.loader.featurizer.structure_featurizer`` -> ``PlinderDataset``).
    The upstream DiffDock conformer-generation, torsion-optimization and
    receptor-graph helpers were removed here as dead code.
"""

allowable_features = {
    "possible_atomic_num_list": list(range(1, 119)) + ["misc"],
    "possible_chirality_list": [
        "CHI_UNSPECIFIED",
        "CHI_TETRAHEDRAL_CW",
        "CHI_TETRAHEDRAL_CCW",
        "CHI_OTHER",
    ],
    "possible_degree_list": [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, "misc"],
    "possible_numring_list": [0, 1, 2, 3, 4, 5, 6, "misc"],
    "possible_implicit_valence_list": [0, 1, 2, 3, 4, 5, 6, "misc"],
    "possible_formal_charge_list": [-5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, "misc"],
    "possible_numH_list": [0, 1, 2, 3, 4, 5, 6, 7, 8, "misc"],
    "possible_number_radical_e_list": [0, 1, 2, 3, 4, "misc"],
    "possible_hybridization_list": ["SP", "SP2", "SP3", "SP3D", "SP3D2", "misc"],
    "possible_is_aromatic_list": [False, True],
    "possible_is_in_ring3_list": [False, True],
    "possible_is_in_ring4_list": [False, True],
    "possible_is_in_ring5_list": [False, True],
    "possible_is_in_ring6_list": [False, True],
    "possible_is_in_ring7_list": [False, True],
    "possible_is_in_ring8_list": [False, True],
}


def safe_index(l, e):
    """Return index of element e in list l. If e is not present, return the last index"""
    try:
        return l.index(e)
    except:
        return len(l) - 1


def lig_atom_featurizer(mol: Chem.rdchem.Mol) -> list[list[int]]:
    ringinfo = mol.GetRingInfo()
    atom_features_list = []
    for idx, atom in enumerate(mol.GetAtoms()):
        chiral_tag = str(atom.GetChiralTag())
        if chiral_tag in [
            "CHI_SQUAREPLANAR",
            "CHI_TRIGONALBIPYRAMIDAL",
            "CHI_OCTAHEDRAL",
        ]:
            chiral_tag = "CHI_OTHER"

        atom_features_list.append(
            [
                safe_index(
                    allowable_features["possible_atomic_num_list"], atom.GetAtomicNum()
                ),
                allowable_features["possible_chirality_list"].index(str(chiral_tag)),
                safe_index(
                    allowable_features["possible_degree_list"], atom.GetTotalDegree()
                ),
                safe_index(
                    allowable_features["possible_formal_charge_list"],
                    atom.GetFormalCharge(),
                ),
                safe_index(
                    allowable_features["possible_implicit_valence_list"],
                    atom.GetImplicitValence(),
                ),
                safe_index(
                    allowable_features["possible_numH_list"], atom.GetTotalNumHs()
                ),
                safe_index(
                    allowable_features["possible_number_radical_e_list"],
                    atom.GetNumRadicalElectrons(),
                ),
                safe_index(
                    allowable_features["possible_hybridization_list"],
                    str(atom.GetHybridization()),
                ),
                allowable_features["possible_is_aromatic_list"].index(
                    atom.GetIsAromatic()
                ),
                safe_index(
                    allowable_features["possible_numring_list"],
                    ringinfo.NumAtomRings(idx),
                ),
                allowable_features["possible_is_in_ring3_list"].index(
                    ringinfo.IsAtomInRingOfSize(idx, 3)
                ),
                allowable_features["possible_is_in_ring4_list"].index(
                    ringinfo.IsAtomInRingOfSize(idx, 4)
                ),
                allowable_features["possible_is_in_ring5_list"].index(
                    ringinfo.IsAtomInRingOfSize(idx, 5)
                ),
                allowable_features["possible_is_in_ring6_list"].index(
                    ringinfo.IsAtomInRingOfSize(idx, 6)
                ),
                allowable_features["possible_is_in_ring7_list"].index(
                    ringinfo.IsAtomInRingOfSize(idx, 7)
                ),
                allowable_features["possible_is_in_ring8_list"].index(
                    ringinfo.IsAtomInRingOfSize(idx, 8)
                ),
            ]
        )
    return atom_features_list
