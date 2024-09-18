# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
# mypy: disable-error-code="no-untyped-def, no-untyped-call, attr-defined, assignment, arg-type, var-annotated"
# ruff: noqa

import copy

import networkx as nx
import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, GetPeriodicTable, rdMolTransforms
from rdkit.Chem.rdchem import BondType as BT
from scipy.optimize import differential_evolution

from plinder.core.utils.log import setup_logger

RDLogger.DisableLog("rdApp.*")
LOG = setup_logger(__name__)

"""
    Preprocessing functions taken or adapted from https://github.com/gcorso/DiffDock
"""


def GetDihedral(conf, atom_idx):
    return rdMolTransforms.GetDihedralRad(
        conf, atom_idx[0], atom_idx[1], atom_idx[2], atom_idx[3]
    )


def SetDihedral(conf, atom_idx, new_vale):
    rdMolTransforms.SetDihedralRad(
        conf, atom_idx[0], atom_idx[1], atom_idx[2], atom_idx[3], new_vale
    )


def apply_changes(mol, values, rotatable_bonds, conf_id):
    opt_mol = copy.copy(mol)
    [
        SetDihedral(opt_mol.GetConformer(conf_id), rotatable_bonds[r], values[r])
        for r in range(len(rotatable_bonds))
    ]
    return opt_mol


def optimize_rotatable_bonds(
    mol,
    true_mol,
    rotatable_bonds,
    probe_id=-1,
    ref_id=-1,
    seed=0,
    popsize=15,
    maxiter=500,
    mutation=(0.5, 1),
    recombination=0.8,
):
    opt = OptimizeConformer(
        mol, true_mol, rotatable_bonds, seed=seed, probe_id=probe_id, ref_id=ref_id
    )
    max_bound = [np.pi] * len(opt.rotatable_bonds)
    min_bound = [-np.pi] * len(opt.rotatable_bonds)
    bounds = (min_bound, max_bound)
    bounds = list(zip(bounds[0], bounds[1]))

    # Optimize conformations
    result = differential_evolution(
        opt.score_conformation,
        bounds,
        maxiter=maxiter,
        popsize=popsize,
        mutation=mutation,
        recombination=recombination,
        disp=False,
        seed=seed,
    )
    opt_mol = apply_changes(opt.mol, result["x"], opt.rotatable_bonds, conf_id=probe_id)

    return opt_mol


class OptimizeConformer:
    def __init__(
        self, mol, true_mol, rotatable_bonds, probe_id=-1, ref_id=-1, seed=None
    ):
        super(OptimizeConformer, self).__init__()
        if seed:
            np.random.seed(seed)
        self.rotatable_bonds = rotatable_bonds
        self.mol = mol
        self.true_mol = true_mol
        self.probe_id = probe_id
        self.ref_id = ref_id

    def score_conformation(self, values):
        for i, r in enumerate(self.rotatable_bonds):
            SetDihedral(self.mol.GetConformer(self.probe_id), r, values[i])
        return AllChem.AlignMol(self.mol, self.true_mol, self.probe_id, self.ref_id)


def get_torsion_angles(mol):
    torsions_list = []
    G = nx.Graph()
    for i, atom in enumerate(mol.GetAtoms()):
        G.add_node(i)
    nodes = set(G.nodes())
    for bond in mol.GetBonds():
        start, end = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        G.add_edge(start, end)
    for e in G.edges():
        G2 = copy.deepcopy(G)
        G2.remove_edge(*e)
        if nx.is_connected(G2):
            continue
        l = list(sorted(nx.connected_components(G2), key=len)[0])
        if len(l) < 2:
            continue
        n0 = list(G2.neighbors(e[0]))
        n1 = list(G2.neighbors(e[1]))
        torsions_list.append((n0[0], e[0], e[1], n1[0]))
    return torsions_list


def read_molecule(
    molecule_file, sanitize=False, calc_charges=False, remove_hs=False
) -> Chem.rdchem.Mol:
    if molecule_file.endswith(".mol2"):
        mol = Chem.MolFromMol2File(molecule_file, sanitize=False, removeHs=False)
    elif molecule_file.endswith(".sdf"):
        supplier = Chem.SDMolSupplier(molecule_file, sanitize=False, removeHs=False)
        mol = supplier[0]
    elif molecule_file.endswith(".pdbqt"):
        with open(molecule_file) as file:
            pdbqt_data = file.readlines()
        pdb_block = ""
        for line in pdbqt_data:
            pdb_block += "{}\n".format(line[:66])
        mol = Chem.MolFromPDBBlock(pdb_block, sanitize=False, removeHs=False)
    elif molecule_file.endswith(".pdb"):
        mol = Chem.MolFromPDBFile(molecule_file, sanitize=False, removeHs=False)
    else:
        raise ValueError(
            "Expect the format of the molecule_file to be "
            "one of .mol2, .sdf, .pdbqt and .pdb, got {}".format(molecule_file)
        )

    try:
        if sanitize or calc_charges:
            Chem.SanitizeMol(mol)

        if calc_charges:
            # Compute Gasteiger charges on the molecule.
            try:
                AllChem.ComputeGasteigerCharges(mol)
            except:
                LOG.warning("Unable to compute charges for the molecule.")

        if remove_hs:
            mol = Chem.RemoveHs(mol, sanitize=sanitize)

    except Exception:
        # Print stacktrace
        import traceback

        msg = traceback.format_exc()
        LOG.warning(f"Failed to process molecule: {molecule_file}\n{msg}")
        return None

    return mol


def generate_conformer(mol):
    ps = AllChem.ETKDGv2()
    failures, id = 0, -1
    while failures < 3 and id == -1:
        if failures > 0:
            LOG.debug(f"rdkit coords could not be generated. trying again {failures}.")
        id = AllChem.EmbedMolecule(mol, ps)
        failures += 1
    if id == -1:
        LOG.info(
            "rdkit coords could not be generated without using random coords. using random coords now."
        )
        ps.useRandomCoords = True
        AllChem.EmbedMolecule(mol, ps)
        AllChem.MMFFOptimizeMolecule(mol, confId=0)
        return True
    return False


periodic_table = GetPeriodicTable()
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
    "possible_amino_acids": [
        "ALA",
        "ARG",
        "ASN",
        "ASP",
        "CYS",
        "GLN",
        "GLU",
        "GLY",
        "HIS",
        "ILE",
        "LEU",
        "LYS",
        "MET",
        "PHE",
        "PRO",
        "SER",
        "THR",
        "TRP",
        "TYR",
        "VAL",
        "HIP",
        "HIE",
        "TPO",
        "HID",
        "LEV",
        "MEU",
        "PTR",
        "GLV",
        "CYT",
        "SEP",
        "HIZ",
        "CYM",
        "GLM",
        "ASQ",
        "TYS",
        "CYX",
        "GLZ",
        "misc",
    ],
    "possible_atom_type_2": [
        "C*",
        "CA",
        "CB",
        "CD",
        "CE",
        "CG",
        "CH",
        "CZ",
        "N*",
        "ND",
        "NE",
        "NH",
        "NZ",
        "O*",
        "OD",
        "OE",
        "OG",
        "OH",
        "OX",
        "S*",
        "SD",
        "SG",
        "misc",
    ],
    "possible_atom_type_3": [
        "C",
        "CA",
        "CB",
        "CD",
        "CD1",
        "CD2",
        "CE",
        "CE1",
        "CE2",
        "CE3",
        "CG",
        "CG1",
        "CG2",
        "CH2",
        "CZ",
        "CZ2",
        "CZ3",
        "N",
        "ND1",
        "ND2",
        "NE",
        "NE1",
        "NE2",
        "NH1",
        "NH2",
        "NZ",
        "O",
        "OD1",
        "OD2",
        "OE1",
        "OE2",
        "OG",
        "OG1",
        "OH",
        "OXT",
        "SD",
        "SG",
        "misc",
    ],
}
bonds = {BT.SINGLE: 0, BT.DOUBLE: 1, BT.TRIPLE: 2, BT.AROMATIC: 3}

lig_feature_dims = (
    list(
        map(
            len,
            [
                allowable_features["possible_atomic_num_list"],
                allowable_features["possible_chirality_list"],
                allowable_features["possible_degree_list"],
                allowable_features["possible_formal_charge_list"],
                allowable_features["possible_implicit_valence_list"],
                allowable_features["possible_numH_list"],
                allowable_features["possible_number_radical_e_list"],
                allowable_features["possible_hybridization_list"],
                allowable_features["possible_is_aromatic_list"],
                allowable_features["possible_numring_list"],
                allowable_features["possible_is_in_ring3_list"],
                allowable_features["possible_is_in_ring4_list"],
                allowable_features["possible_is_in_ring5_list"],
                allowable_features["possible_is_in_ring6_list"],
                allowable_features["possible_is_in_ring7_list"],
                allowable_features["possible_is_in_ring8_list"],
            ],
        )
    ),
    0,
)  # number of scalar features

rec_atom_feature_dims = (
    list(
        map(
            len,
            [
                allowable_features["possible_amino_acids"],
                allowable_features["possible_atomic_num_list"],
                allowable_features["possible_atom_type_2"],
                allowable_features["possible_atom_type_3"],
            ],
        )
    ),
    0,
)

rec_residue_feature_dims = (
    list(map(len, [allowable_features["possible_amino_acids"]])),
    0,
)


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
                # g_charge if not np.isnan(g_charge) and not np.isinf(g_charge) else 0.
            ]
        )
    return atom_features_list


def get_lig_graph(mol):
    atom_feats = lig_atom_featurizer(mol)

    row, col, edge_type = [], [], []
    for bond in mol.GetBonds():
        start, end = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        row += [start, end]
        col += [end, start]
        edge_type += (
            2 * [bonds[bond.GetBondType()]]
            if bond.GetBondType() != BT.UNSPECIFIED
            else [0, 0]
        )
    if mol.GetNumConformers() > 0:
        lig_coords = mol.GetConformer().GetPositions()
    return


def get_lig_graph_with_matching(
    mol_,
    popsize=15,
    maxiter=15,
    num_conformers=1,
    remove_hs=False,
    tries=10,
    skip_matching=False,
) -> None:
    mol_maybe_noh = copy.deepcopy(mol_)
    if remove_hs:
        mol_maybe_noh = Chem.RemoveHs(mol_maybe_noh, sanitize=True)
        mol_maybe_noh = AllChem.RemoveAllHs(mol_maybe_noh)

    _tmp = copy.deepcopy(mol_)
    if remove_hs:
        _tmp = Chem.RemoveHs(_tmp, sanitize=True)
    _tmp = AllChem.RemoveAllHs(_tmp)
    rotatable_bonds = get_torsion_angles(_tmp)

    for i in range(num_conformers):
        mols, rmsds = [], []
        for _ in range(tries):
            mol_rdkit = copy.deepcopy(mol_)

            mol_rdkit.RemoveAllConformers()
            mol_rdkit = AllChem.AddHs(mol_rdkit)
            generate_conformer(mol_rdkit)
            if remove_hs:
                mol_rdkit = Chem.RemoveHs(mol_rdkit, sanitize=True)
            mol_rdkit = AllChem.RemoveAllHs(mol_rdkit)
            mol = AllChem.RemoveAllHs(copy.deepcopy(mol_maybe_noh))
            if rotatable_bonds and not skip_matching:
                optimize_rotatable_bonds(
                    mol_rdkit,
                    mol,
                    rotatable_bonds,
                    popsize=popsize,
                    maxiter=maxiter,
                )
            mol.AddConformer(mol_rdkit.GetConformer())
            rms_list = []
            AllChem.AlignMolConformers(mol, RMSlist=rms_list)
            mol_rdkit.RemoveAllConformers()
            mol_rdkit.AddConformer(mol.GetConformers()[1])
            mols.append(mol_rdkit)
            rmsds.append(rms_list[0])
        mol_rdkit = mols[np.argmin(rmsds)]
        if i == 0:
            rmsd_matching = min(rmsds)
            get_lig_graph(mol_rdkit)
        else:
            pos = mol_rdkit.GetConformer().GetPositions()
    return
