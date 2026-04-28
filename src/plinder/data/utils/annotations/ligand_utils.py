# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
from __future__ import annotations

import itertools
import logging
import re
import sqlite3
import typing as ty
from collections import Counter, defaultdict
from functools import cache, cached_property
from pathlib import Path

import biotite.structure as struc
import biotite.structure.info as bt_info
import numpy as np
import pandas as pd
from pydantic import BeforeValidator, Field
from rdkit import Chem, RDLogger
from rdkit.Chem import QED, Crippen, rdMolDescriptors
from rdkit.Chem import rdMolDescriptors as rdMD
from rdkit.Chem.rdchem import Mol

from plinder.core.utils.config import get_config
from plinder.core.utils.constants import BASE_DIR
from plinder.data.utils.annotations.interaction_utils import (
    extract_ligand_links_to_neighbouring_chains,
    run_peppr_interactions,
)
from plinder.data.utils.annotations.interface_gap import (
    annotate_interface_gaps_per_chain,
)
from plinder.data.utils.annotations.protein_utils import Chain, sequences_match_core
from plinder.data.utils.annotations.utils import DocBaseModel

_PRD_DB_PATH = str(BASE_DIR / "data/utils/annotations/static_files/prdcc.chemlib")
LOG = logging.getLogger(__name__)


def _template_from_user_smiles(
    comp_id: str,
    smiles: str,
    cif_atom_names: list[str],
) -> "Chem.Mol | None":
    """Build a stereo-assigned template Mol from a user-supplied SMILES.

    Used as a CCD fallback when a custom residue (e.g. Boltz ``LIG``) is
    not in the Chemical Component Dictionary. Assumes the SMILES
    heavy-atom parse order matches the CIF heavy-atom order — the same
    positional convention used by :func:`assign_bond_orders_from_smiles`.

    Stereo is assigned from SMILES parity tags (``@``/``@@``) directly,
    no 3D embed needed. PDB atom names from the CIF are stamped onto the
    template atoms so :func:`compare_stereo_to_template` can match by name.

    Returns ``None`` if the SMILES can't be parsed or the heavy-atom
    count disagrees with the CIF (the caller then falls back to ``None``
    for stereo_matches, matching pre-existing behaviour).
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    mol = Chem.RemoveHs(mol, sanitize=False)
    if mol.GetNumAtoms() != len(cif_atom_names):
        LOG.warning(
            f"_template_from_user_smiles: atom count mismatch for {comp_id} "
            f"({mol.GetNumAtoms()} in SMILES vs {len(cif_atom_names)} in CIF) — "
            "skipping SMILES-based stereo check"
        )
        return None
    Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
    for atom, atom_name in zip(mol.GetAtoms(), cif_atom_names):
        info = Chem.AtomPDBResidueInfo()
        info.SetName(atom_name)
        info.SetResidueName(comp_id)
        info.SetResidueNumber(1)
        atom.SetMonomerInfo(info)
    return mol


def _check_stereo_vs_template(
    resolved_mol: "Chem.Mol",
    custom_templates: dict[str, "Chem.Mol"] | None = None,
) -> bool | None:
    """Compare resolved 3D stereo against a stereo template per residue.

    Template source precedence:
      1. ``custom_templates[resname]`` if provided — user-supplied SMILES
         templates win over CCD because the caller explicitly knows CCD
         is wrong or missing (biotite ships a generic placeholder for
         some codes like ``LIG`` that would otherwise silently hide
         stereo mismatches).
      2. :func:`_get_ccd_mol(resname)` — CCD ideal coordinates.
      3. Return ``None`` for this residue if neither source yields a
         template.

    Delegates to :func:`compare_stereo_to_template` for the actual CIP
    comparison. Handles multi-residue ligands (e.g. glycans) by
    checking each residue copy independently.

    Returns ``True`` if all residues match or are achiral, ``False`` if
    any stereo mismatch, ``None`` if no template was available for any
    residue.
    """
    from plinder.core.structure.smallmols_utils import compare_stereo_to_template

    # Group atoms by (resname, res_id) to handle repeated residue names
    residue_atoms: dict[tuple[str, int], list[int]] = {}
    for atom in resolved_mol.GetAtoms():
        info = atom.GetPDBResidueInfo()
        if info is None:
            raise ValueError(
                f"Atom {atom.GetIdx()} in resolved mol has no PDB residue info"
            )
        key = (info.GetResidueName().strip(), info.GetResidueNumber())
        residue_atoms.setdefault(key, []).append(atom.GetIdx())

    results: list[bool | None] = []
    for (resname, res_id), atom_indices in residue_atoms.items():
        # User-supplied custom templates take precedence over CCD: if
        # the caller provided a SMILES template for this residue, they
        # explicitly know CCD is wrong or missing (biotite ships a
        # generic placeholder for some codes like "LIG" that would
        # otherwise hide stereo mismatches).
        template_mol = None
        if custom_templates is not None:
            template_mol = custom_templates.get(resname)
        if template_mol is None:
            template_mol = _get_ccd_mol(resname)
        if template_mol is None:
            results.append(None)
            continue

        frag = Chem.RWMol(resolved_mol)
        remove = [
            a.GetIdx()
            for a in resolved_mol.GetAtoms()
            if a.GetIdx() not in atom_indices
        ]
        frag.BeginBatchEdit()
        for idx in sorted(remove, reverse=True):
            frag.RemoveAtom(idx)
        frag.CommitBatchEdit()

        try:
            results.append(compare_stereo_to_template(frag.GetMol(), template_mol))
        except Exception as e:
            LOG.warning(f"Stereo comparison failed for {resname}:{res_id}: {e}")
            results.append(None)

    if not results:
        return None
    if any(r is False for r in results):
        return False
    if any(r is True for r in results):
        return True
    return None


@cache
def _get_ccd_mol(comp_id: str) -> "Chem.Mol | None":
    """Return an RDKit Mol from CCD ideal coordinates with stereo assigned."""
    from plinder.data.utils.annotations.cif_utils import atoms_to_rdkit_mol

    try:
        # biotite has no chiral tags, so atoms_to_rdkit_mol uses 3D to
        # assign stereo via AssignStereochemistryFrom3D
        return atoms_to_rdkit_mol(bt_info.residue(comp_id))
    except Exception as e:
        LOG.warning(f"Failed to get CCD mol for {comp_id}: {e}")
        return None


def _get_ccd_smiles(comp_id: str) -> str | None:
    """Get SMILES from CCD via biotite, with stereochemistry from ideal 3D."""
    mol = _get_ccd_mol(comp_id)
    if mol is None:
        return None
    return str(Chem.MolToSmiles(mol))


def _get_prd_smiles(comp_id: str) -> str | None:
    """Get SMILES from PRD library (SQLite)."""
    try:
        conn = sqlite3.connect(_PRD_DB_PATH)
        cursor = conn.cursor()
        cursor.execute("SELECT smiles FROM chem_compounds WHERE tlc = ?", (comp_id,))
        row = cursor.fetchone()
        conn.close()
        if row and row[0]:
            return str(row[0])
    except Exception:
        LOG.warning(f"Failed to fetch PRD SMILES for {comp_id}")
    return None


def lig_has_dummies(
    ligand_code: str,
    dummy_lig_list: list[str] = [
        "DUM",
        "UNX",
        "ASX",
        "GLX",
        "UNL",
        "UNK",
        "UPL",
        "DN",
        "N",
    ],
) -> bool:
    """Check for ccd codes containing dummy/unknown entries

    Args:
        ligand_code str: ligand CCD code
        dummy_lig_list (list, optional): list of ccd codes for unknown or dummy entries.
        Defaults to ['DUM', 'UNX', 'UNL', 'UNK', 'UPL', 'DN', 'N'].

    Returns:
        bool: if ligand considered as dummy and treated as artifact
    """
    # check for dummy list including composites, too!
    return len(set(ligand_code.split("-")).intersection(dummy_lig_list)) > 0


def get_ccd_smiles_dict(ciffile: Path) -> dict[str, str]:
    """Load CCD component SMILES from a parquet file next to *ciffile*."""
    df = pd.read_parquet(ciffile.parent / "components.parquet")
    return dict(zip(df["binder_id"], df["canonical_smiles"]))


def sort_ccd_codes(code_list: list[str]) -> list[str]:
    """Pick long first, then alphabetical letters followed by numbers
    Args:
        code_list (Set[str]): set of CCD strings

    Returns:
        List[str]: list of sorted CCD string set
    """
    code_list = sorted(sorted(code_list), key=len, reverse=True)
    final_list = [code for code in code_list if not re.findall("([0-9])", code[0])] + [
        code for code in code_list if re.findall("([0-9])", code[0])
    ]
    return final_list


@cache
def get_ccd_synonyms(data_dir: Path) -> tuple[list[set[str]], dict[str, str]]:
    """Get Synonym dictonary for CCD SMILES
    CCD smiles dict from download_components_cif
    and get_ccd_smiles_dict

    Returns:
        Dict[str, str]: dictonary mapping synonymous CCD code to preferred one
    """
    from plinder.data.pipeline.io import download_components_cif

    ccd_lib_cifpath = download_components_cif(data_dir=data_dir)
    smidict = get_ccd_smiles_dict(ccd_lib_cifpath)
    ccd_df = pd.DataFrame.from_dict(smidict, orient="index").reset_index()
    # note: SMILES assumed to be CANONICALIZED by OE read by get_ccd_smiles_dict()
    ccd_df.columns = ["ccd_code", "SMILES"]
    # remove dummies
    ccd_df = ccd_df[~ccd_df["ccd_code"].apply(lig_has_dummies)]
    ccd_sets = ccd_df.groupby("SMILES").aggregate(list).reset_index()
    # ccd_sets
    ccd_sets["ccd_synonym_count"] = ccd_sets["ccd_code"].apply(lambda x: len(x))
    ccd_dups = ccd_sets[ccd_sets.ccd_synonym_count > 1].copy()
    list_of_synonym_sets = [set(x) for x in ccd_dups["ccd_code"].to_list()]
    # keep unique_ccd as sorted first entry
    ccd_dups["unique_ccd"] = ccd_dups["ccd_code"].apply(lambda x: sort_ccd_codes(x)[0])
    ccd_dups_exp = ccd_dups.explode("ccd_code")
    ccd_dups_exp.index = ccd_dups_exp["ccd_code"]
    ccd_synonym_dict = ccd_dups_exp["unique_ccd"].to_dict()
    return (list_of_synonym_sets, ccd_synonym_dict)


# lazy evaluate data fetches referenced as module globals
# TODO : clean this up and deduplicate extras with pipeline.io
COFACTORS = None
LIST_OF_CCD_SYNONYMS = None
CCD_SYNONYMS_DICT = None
# instantiate artifact list once and reuse variable
ARTIFACTS = None
KINASE_INHIBITORS = None
BINDING_AFFINITY = None


def add_missed_synonyms(current_set: set[str]) -> set[str]:
    """Expand a set of CCD codes with any known synonyms."""
    assert LIST_OF_CCD_SYNONYMS is not None
    missed_synonyms = [
        x.difference(current_set)
        for x in LIST_OF_CCD_SYNONYMS
        if len(x.intersection(current_set)) > 0
    ]
    return set(itertools.chain(*missed_synonyms)).union(current_set)


def get_unique_ccd_longname(longname: str) -> str:
    """Map a composite CCD code to its canonical synonym form."""
    assert CCD_SYNONYMS_DICT is not None

    if longname.startswith("PRD_"):
        return longname
    else:
        return "-".join([CCD_SYNONYMS_DICT.get(s, s) for s in longname.split("-")])


def get_chain_type(chain_type_str: str) -> str:
    """Classify chain type string into ligand category."""
    ct = chain_type_str.lower()
    if "non-polymer" in ct:
        return "SMALLMOLECULE"
    if "polypeptide" in ct:
        return "PEPTIDE"
    if "polydeoxyribonucleotide" in ct and "polyribonucleotide" in ct:
        return "MIXED"
    if "polydeoxyribonucleotide" in ct:
        return "DNA"
    if "polyribonucleotide" in ct:
        return "RNA"
    if "polysaccharide" in ct or "oligosaccharide" in ct or "branched" in ct:
        return "SACCHARIDE"
    if "macrolide" in ct or "cyclic-pseudo-peptide" in ct:
        return "MACROCYCLES"
    return "UNKNOWN"


@cache
def parse_cofactors(data_dir: Path) -> set[str]:
    """Download and parse cofactors.

    Returns
    -------
    Set[str]
        Set of cofactors

    """
    from plinder.data.pipeline.io import download_cofactors

    cofactors_json = download_cofactors(data_dir=data_dir)
    extra = {
        "Ascorbic acid": ["UU3"],
        "Coenzyme F420": ["6J4", "F42"],
        "Factor F430": ["F43", "M43"],
        "Pantetheine": ["PNY"],
        "Pantothenic acids": ["66S", "8Q1", "PAU"],
        "Nicotinamide": ["NCA"],
        "Adenosine nucleotides": [
            "A",
            "AMP",
            "ATP",
            "ADP",
        ],  # + ["ANP"],  # ANP is a mimic-inhibitor
        "Guanosine nucleotides": [
            "G",
            "GTP",
            "GDP",
            "GMP",
            "CPG",
            "G25",
            "5GP",
        ],  # + ["GNP", "GTN"],  # GNP/GTN is inhibitor
        "Cytidine nucleotides": ["C", "C5P", "C25", "CDP", "CTP"],
        "Thymidine nucleotides": ["T", "TMP", "DT", "TTP", "THM", "TYD"],
        "Uridine nucleotides": ["U", "DU", "U5P", "U25", "UMP", "UDP", "UTP"],
        "MIO": ["CRW"],
        "NAD": ["NAH"],
        "Glutathione": ["CYP"],
        "Biopterin": ["HBL", "BH4", "THB"],
        "Tetrahydrofolic acid": ["MEF"],
        "Lumazine": ["DLZ"],
        "Menaquinone": ["MQ8", "MQ9", "MQE", "7MQ"],
        "Heme": ["1CP", "CP3", "MMP", "UP2", "UP3"],
        "Methanopterin": ["H4M", "H4Z"],
        "Lipoamide": ["LPM"],
        "Ubiquinone": ["DCQ", "HQE", "PLQ"],
        "Pyridoxal": ["PXL", "UEG"],
        "Siderophores": ["488", "EB4", "SE8"],
        "Methanofuran": ["MFN"],
        "Vitamin A": ["BCR", "ECH", "EQ3", "RAW"],
        "Vitamin K1": ["PQN"],
        "CHLOROPHYLL and similar": [
            "CLA",
            "CHL",
            "CL0",
            "CL1",
            "CL2",
            "CL7",
            "BCB",
            "BCL",
            "07D",
            "G9R",
            "PEB",
            "PUB",
            "CYC",
            "BPH",
        ],
        # "Lipids": ["SPH"], # TODO: ?
        # "Sugars": ["NAG", "BCG", "GLC"], # TODO: more?
        #
    }
    cofactors = set()
    for c in cofactors_json:
        for c_list in cofactors_json[c]:
            cofactors |= set(c_list.get("cofactors", []))
    for c in extra:
        cofactors |= set(extra[c])

    # add missed synonyms
    cofactors = add_missed_synonyms(cofactors)

    return cofactors


@cache
def parse_artifacts() -> set[str]:
    """Get and parse artifacts
    Returns:
        set[str]: set[str]
    """
    artifact_log = BASE_DIR / "utils/annotations/static_files/artifacts_badlist.csv"
    with open(artifact_log, "r") as f:
        lines = f.readlines()
    artifacts = {l.strip() for l in lines if not l.startswith("#")}
    # add missed synonyms
    artifacts = add_missed_synonyms(artifacts)
    return artifacts


@cache
def parse_kinase_inhibitors(data_dir: Path) -> set[str]:
    """Load set of CCD codes for known kinase inhibitors."""
    from plinder.data.pipeline.io import download_kinase_data

    kinase_ligand_path = download_kinase_data(data_dir=data_dir)
    kinase_ligand_path = kinase_ligand_path.with_name("kinase_ligand_ccd_codes.parquet")
    kinase_ligand_df = pd.read_parquet(kinase_ligand_path)
    return set(kinase_ligand_df["PDB-code"])


@cache
def get_binding_affinity(data_dir: Path) -> ty.Any:
    """Load BindingDB affinity data (pchembl values + target sequences)."""
    from plinder.data.pipeline.io import download_affinity_data

    return download_affinity_data(data_dir=data_dir)


def get_num_resolved_heavy_atoms(resolved_smiles: str) -> int:
    """Count heavy atoms in the resolved SMILES (0 if unparseable)."""
    matched_mol = Chem.MolFromSmiles(resolved_smiles, sanitize=False)
    if matched_mol is None:
        return 0
    return int(rdMD.CalcNumHeavyAtoms(matched_mol))


def get_len_of_longest_linear_hydrocarbon_linker(
    mol: Mol,
    max_count: int = 50,
    link_unit_smarts: str = "[#6D2R0]",
) -> int:
    """Estimate maximum linker length defined by link_unit_smarts, eg.
    unbranched hydrocarbons (default)

    Args:
        mol (Mol): RDKit molecule
        max_count (int, optional):
            Max count for linker. Defaults to 50.
        link_unit_smarts (str, optional):
            Linker unit defined by SMARTS. Defaults to "[#6D2R0]".

    Returns:
        int: maximum linker length defined by link_unit_smarts (default: unbranched hydrocarbon)
    """
    try:
        # needs ring info!
        Chem.SanitizeMol(
            mol, sanitizeOps=Chem.rdmolops.SanitizeFlags.SANITIZE_SYMMRINGS
        )
        # length of longest hydrocarbon chain (excludes the ends and rings)
        for i in range(max_count):
            # chain_smarts = "[#6D2R0,#6D1R0]" * (i+1) # includes the ends
            chain_smarts = "~".join([link_unit_smarts] * (i + 1))
            if len(mol.GetSubstructMatches(Chem.MolFromSmarts(chain_smarts))) == 0:
                return i
        # TODO: what to do if fails or not found? now returns -1
        return max_count + 100
    except Exception as e:
        logging.warning(
            f"Error in calculating longest linear hydrocarbon linker for {mol.GetProp('_Name')}: {e}"
        )
        return max_count + 100


def is_excluded_mol(
    smiles: str,
    min_C_threshold: int = 2,
    min_HA_threshold: int = 5,
    max_charge: int = 2,
    max_linear_hydrocarbon_linker: int = 12,
) -> bool:
    """Exclude some molecules by default as useless for druglikeness
    Uses OR logic for violating rules:
        - less than 2 carbon atoms
        - less than 5 non-hydrogen atoms
        - charge larger than +/- 2
        - unbranched hydrocarbon linker no longer than 12

    Args:
        smiles (str): molecule SMILES
        min_C_threshold (int, optional):
            Minimum carbon atom count. Defaults to 2.
        min_HA_threshold (int, optional):
            Minimum non-hydrogen atom count. Defaults to 5.
        max_charge (int, optional):
            Maximum allowed absolute charge. Defaults to 2.
        max_linear_hydrocarbon_linker (int, optional):
            Maximum allowed unbranched hydrocarbon linker. Defaults to 12.

    Returns:
        bool: should molecule be considered as artifact
    """
    mol = Chem.MolFromSmiles(smiles, sanitize=False)

    # get heavy atom and carbon counts
    carbon = Chem.MolFromSmarts("[#6]")
    numC = len(mol.GetSubstructMatches(carbon))
    numHA = mol.GetNumHeavyAtoms()

    if numHA < min_HA_threshold or numC < min_C_threshold:
        return True

    # get formal charge
    charge = Chem.rdmolops.GetFormalCharge(mol)
    if abs(charge) > max_charge:
        return True
    elif (
        get_len_of_longest_linear_hydrocarbon_linker(mol)
        > max_linear_hydrocarbon_linker
    ):
        return True
    else:
        return False


def is_single_atom_or_ion(mol: Mol) -> bool:
    """True if the molecule is a single non-organic heavy atom (metal ion)."""
    numHA = mol.GetNumHeavyAtoms()
    skip_single_elems = Chem.MolFromSmarts("[#6,#1,#0,#7,#8,#15,#16,#34,#52]")
    numCHNOPSetc = len(mol.GetSubstructMatches(skip_single_elems))
    return numHA == 1 and numCHNOPSetc == 0


def validate_chain_residue(obj: dict[str, ty.Any]) -> dict[str, ty.Any]:
    """Recursively coerce string dict keys to ints or tuples for pydantic."""
    clean = {}
    for k, v in obj.items():
        if isinstance(k, str):
            if "," in k:
                key: ty.Any = tuple(k.split(","))
            else:
                try:
                    key = int(k)
                except ValueError:
                    key = k
        else:
            key = k
        if isinstance(v, dict):
            clean[key] = validate_chain_residue(v)
        else:
            clean[key] = v
    return clean


CrystalContacts = ty.Annotated[
    dict[tuple[str, int], set[int]],
    BeforeValidator(validate_chain_residue),
    Field(default_factory=dict),
]


class Ligand(DocBaseModel):
    pdb_id: str = Field(
        default_factory=str,
        description="__RCSB PDB ID, see https://mmcif.wwpdb.org/dictionaries/mmcif_pdbx_v50.dic/Items/_entry.id.html",
    )
    biounit_id: str = Field(default_factory=str, description="__Biounit id")
    asym_id: str = Field(default_factory=str, description="Ligand chain asymmetric id")
    instance: int = Field(default_factory=int, description="Biounit instance ID")
    ccd_code: str = Field(
        default_factory=str,
        description="Ligand Chemical Component Dictionary (CCD) code",
    )
    # TODO: rename plip_type → chain_type; name kept for backward compatibility
    # (PLIP tool is no longer used — replaced by peppr)
    plip_type: str = Field(
        default_factory=str, description="Ligand chain type classification"
    )
    bird_id: str = Field(default_factory=str, description="Ligand BIRD id")
    centroid: list[float] = Field(
        default_factory=list, description="Ligand center of geometry"
    )
    smiles: str = Field(
        default_factory=str,
        description="Ligand SMILES from CCD/PRD lookup, or derived from resolved 3D if not in dictionary",
    )
    resolved_smiles: str = Field(
        default_factory=str,
        description="SMILES from resolved 3D coordinates: bond orders from CCD template, stereochemistry from 3D geometry",
    )
    resolved_stereo_matches_template: bool | None = Field(
        default=None,
        description="Whether resolved 3D stereo matches CCD template (True if achiral; None if no template)",
    )
    residue_numbers: list[int] = Field(
        default_factory=list, description="__Ligand residue numbers"
    )
    rdkit_canonical_smiles: str | None = Field(
        default=None,
        description="RDKit canonical SMILES (same as smiles; kept for schema compatibility)",
    )
    molecular_weight: float | None = Field(default=None, description="Molecular weight")
    crippen_clogp: float | None = Field(
        default=None,
        description="Ligand Crippen MlogP, see https://www.rdkit.org/docs/source/rdkit.Chem.Crippen.html",
    )
    num_rot_bonds: int | None = Field(
        default=None, description="Number of rotatable bonds"
    )
    num_hbd: int | None = Field(
        default=None, description="Number of hydrogen bond donors"
    )
    num_hba: int | None = Field(
        default=None, description="Number of hydrogen bond acceptors"
    )
    num_rings: int | None = Field(default=None, description="Number of rings")
    num_heavy_atoms: int | None = Field(
        default=None, description="Number of heavy atoms"
    )
    is_covalent: bool = Field(
        default=False, description="Indicator of whether a ligand  is a covalent ligand"
    )
    covalent_linkages: set[str] = Field(
        default_factory=set[str],
        description="Ligand covalent linkages from _struct_conn (conn_type_id='covale'), "
        + "format: {auth_seq}:{comp_id}:{chain}:{seq}:{atom}__{auth_seq}:{comp_id}:{chain}:{seq}:{atom}",
    )
    neighboring_residues: dict[str, list[int]] = Field(
        default_factory=dict,
        description="Dictionary of neighboring residues, with {instance}.{chain} key and residue number value",
    )
    neighboring_ligands: list[str] = Field(
        default_factory=list,
        description="__List of neighboring ligands {instance}.{chain}",
    )
    receptor_seqres: dict[str, str] = Field(
        default_factory=dict,
        description="__SEQRES sequences of neighboring receptor chains for affinity validation",
    )
    interacting_residues: dict[str, list[int]] = Field(
        default_factory=dict,
        description="Dictionary of interacting residues, with {instance}.{chain} key and residue number value",
    )
    interacting_ligands: list[str] = Field(
        default_factory=list,
        description="__List of interacting ligands {instance}.{chain}",
    )
    # TODO: rename interactions description; hash format kept for backward compatibility
    # (now computed by peppr, not PLIP)
    interactions: dict[str, dict[int, list[str]]] = Field(
        default_factory=dict,
        description="__Dictionary of {instance}.{chain} to residue number to list of interaction hashes",
    )
    neighboring_residue_threshold: float = Field(
        default=6.0,
        description="__Maximum distance to consider receptor residues (protein/NA) neighboring",
    )
    neighboring_ligand_threshold: float = Field(
        default=4.0, description="__Maximum distance to consider ligands neighboring"
    )
    num_neighboring_ppi_atoms_within_4A_of_gap: int | None = Field(
        default=None,
        description="Number of missing neighboring protein-protein interface atoms within 4 Å of ligand of interest",
    )
    num_neighboring_ppi_atoms_within_8A_of_gap: int | None = Field(
        default=None,
        description="Number of missing neighboring protein-protein interface atoms within 8 Å of ligand of interest",
    )
    num_missing_ppi_interface_residues: int | None = Field(
        default=None,
        description="Number of missing neighboring protein-protein interface residues within 4 Å of ligand of interest",
    )
    num_pli_atoms_within_4A_of_gap: int | None = Field(
        default=None,
        description="Number of missing neighboring protein-ligand interface atoms within 4 Å of ligand of interest",
    )
    num_pli_atoms_within_8A_of_gap: int | None = Field(
        default=None,
        description="Number of missing neighboring protein-ligand interface atoms within 8 Å of ligand of interest",
    )
    num_missing_pli_interface_residues: int | None = Field(
        default=None,
        description="Number of missing neighboring protein-ligand interface residues within 4 Å of ligand of interest",
    )
    num_resolved_heavy_atoms: int | None = Field(
        default=None, description="Number of resolved heavy atoms in a ligand"
    )
    num_unresolved_heavy_atoms: int | None = Field(
        default=None, description="Number of unresolved heavy atoms in a ligand"
    )
    tpsa: float | None = Field(
        default=None, description="Topological polar surface area"
    )
    qed: float | None = Field(
        default=None,
        description="Ligand QED score, a measure of drug-likeness, see https://www.rdkit.org/new_docs/source/rdkit.Chem.QED.html",
    )
    is_ion: bool = Field(
        default=False, description="Indicator of whether a ligand  is an ion"
    )
    is_lipinski: bool = Field(
        default=False,
        description="Indicator of whether a ligand satisfies Lipinski Ro5",
    )
    is_fragment: bool = Field(
        default=False,
        description="Indicator of whether a ligand satisfies fragment Ro3",
    )
    is_oligo: bool = Field(
        default=False,
        description="Indicator of whether a ligand  is an oligopeptide, oligosaccharide or oligopeptide",
    )
    is_cofactor: bool = Field(
        default=False, description="Indicator of whether a ligand is a cofactor"
    )
    in_artifact_list: bool = Field(
        default=False,
        description="Indicator of whether a ligand is in the artifact list",
    )
    is_artifact: bool = Field(
        default=False, description="Indicator of whether a ligand is an artifact"
    )
    is_other: bool = Field(
        default=False,
        description="Indicator of whether a ligand type is not classified as any types of small molecule "
        + "(Lipinski, Fragment or covalent), ion, cofactor, oligo (peptide, saccharide or nucleotide) or artifact",
    )
    is_invalid: bool = Field(
        default=False, description="Indicator of whether a ligand is invalid"
    )
    unique_ccd_code: str | None = Field(
        default=None, description="Ligand representative CCD code after de-duplicating"
    )
    crystal_contacts: CrystalContacts = Field(
        default_factory=dict,
        description="__Dictionary of {chain} to residue number to set of interacting crystal contacts",
    )
    waters: dict[str, list[int]] = Field(
        default_factory=dict,
        description="__Dictionary of {instance}.{chain} to list of interacting water residue numbers",
    )
    posebusters_result: dict[str, ty.Any] = Field(
        default_factory=dict,
        description="__Results from running posebusters with 're-dock'",
    )
    """Ligand annotation dataclass.

    Holds structural, chemical, and interaction annotations for a single
    ligand chain in a protein–ligand (or NA–ligand) complex.
    """

    def set_rdkit(self) -> None:
        """Compute RDKit molecular descriptors from ``self.smiles``."""
        try:
            rdkit_compatible_mol = Chem.MolFromSmiles(self.smiles)
            # smiles is already canonical (from MolToSmiles); kept for schema compat
            self.rdkit_canonical_smiles = self.smiles
            self.molecular_weight = rdMolDescriptors.CalcExactMolWt(
                rdkit_compatible_mol
            )
            self.num_rot_bonds = rdMolDescriptors.CalcNumRotatableBonds(
                rdkit_compatible_mol
            )
            self.num_hba = rdMolDescriptors.CalcNumHBA(rdkit_compatible_mol)
            self.num_hbd = rdMolDescriptors.CalcNumHBD(rdkit_compatible_mol)
            self.crippen_clogp = Crippen.MolLogP(rdkit_compatible_mol)
            self.num_rings = rdMolDescriptors.CalcNumRings(rdkit_compatible_mol)
            self.num_heavy_atoms = rdMolDescriptors.CalcNumHeavyAtoms(
                rdkit_compatible_mol
            )
            self.tpsa = rdMolDescriptors.CalcTPSA(rdkit_compatible_mol)
            self.qed = QED.qed(rdkit_compatible_mol)
            self.num_resolved_heavy_atoms = get_num_resolved_heavy_atoms(
                self.resolved_smiles
            )

            if self.num_heavy_atoms and self.num_resolved_heavy_atoms:
                self.num_unresolved_heavy_atoms = (
                    self.num_heavy_atoms - self.num_resolved_heavy_atoms
                )
            # classify ligand based on above molecule
            self.classify_ligand_type(rdkit_compatible_mol)

        except Exception as e:
            logging.warning(f"Error in setting rdkit for {self.id}: {e}")
            # Multi-residue ligands (peptides) may fail SMILES derivation
            # but are still structurally valid
            if self.smiles is None:
                self.is_invalid = True

    def classify_ligand_type(self, mol: Mol) -> None:
        """Classify ligand as ion, Lipinski, fragment, oligo, or artifact.

        Uses SMARTS patterns and Lipinski rules to assign granular type
        beyond the chain-type classification.

        Note
        ----
        Oligo smarts obtained from https://doi.org/10.1021/acs.jcim.3c01573

        Args:
            mol (Mol): RDKit compatible molecule
        """
        RDLogger.DisableLog("rdApp.*")
        oligo_smarts = {
            "oligopeptide": Chem.MolFromSmarts("C(=O)C[N;D2,D3]C(=[O;D1])CN"),
            "oligosaccharide": Chem.MolFromSmarts(
                "O-[C;R0,R1]-[C;R0,R1]-[O,S;R0;D2]-[C;R1]-[O;R1]"
            ),
            "oligonucleotide": Chem.MolFromSmarts("P(=O)([O-,OH])(OC[C;r5])O[C;r5]"),
        }

        if mol is not None:
            if is_single_atom_or_ion(mol):
                self.is_ion = True
            else:
                for sm in oligo_smarts:
                    try:
                        if mol.HasSubstructMatch(oligo_smarts[sm]):
                            # TODO: review - reduced to one class
                            self.is_oligo = True
                    except RuntimeError:
                        self.is_invalid = True
        else:
            self.is_invalid = True

        # Lipinski like Ro3 and Ro5 - for non ions only
        if (
            self.is_ion == False
            and self.is_invalid == False
            and self.molecular_weight is not None
            and self.crippen_clogp is not None
            and self.num_hbd is not None
            and self.num_hba is not None
        ):
            if (
                self.molecular_weight < 300
                and self.crippen_clogp < 3
                and self.num_hbd <= 3
                and self.num_hba <= 3
            ):
                self.is_fragment = True
                self.is_lipinski = True
            elif (
                self.molecular_weight < 500
                and self.crippen_clogp < 5
                and self.num_hbd <= 5
                and self.num_hba <= 10
            ):
                self.is_lipinski = True

    @classmethod
    def from_pli(
        cls,
        pdb_id: str,
        biounit_id: str,
        biounit: ty.Any,
        ligand_instance: int,
        ligand_chain: Chain,
        residue_numbers: list[int],
        ligand_like_chains: dict[str, str],
        interface_proximal_gaps: dict[str, dict[tuple[str, str], dict[str, int]]],
        all_covalent_dict: dict[str, list[tuple[str, str]]],
        # TODO: rename plip_complex_threshold -> complex_threshold
        plip_complex_threshold: float = 10.0,
        neighboring_residue_threshold: float = 6.0,
        neighboring_ligand_threshold: float = 4.0,
        data_dir: ty.Optional[Path] = None,
        chain_to_seqres: dict[str, str] | None = None,
        ligand_smiles_dict: dict[str, str] | None = None,
    ) -> Ligand | None:
        """Build a Ligand from a biounit AtomArray and chain metadata.

        Extracts SMILES (CCD template → resolved 3D fallback), computes
        interactions via peppr, finds neighboring residues/ligands, and
        validates stereochemistry against the CCD template.

        Parameters
        ----------
        pdb_id : str
            PDB entry identifier.
        biounit_id : str
            Biological assembly identifier.
        biounit : struc.AtomArray
            Full biounit atoms with bonds.
        ligand_instance : int
            Instance index within the biounit.
        ligand_chain : Chain
            Chain metadata for the ligand.
        residue_numbers : list[int]
            Residue numbers belonging to this ligand.
        ligand_like_chains : dict[str, str]
            Other ligand-like chains in the entry ``{chain_id: chain_type}``.
        interface_proximal_gaps : dict
            Gap annotation from ``annotate_interface_gaps()``.
        all_covalent_dict : dict[str, list[tuple[str, str]]]
            Covalent linkages by type (``"covale"``, ``"metalc"``, ``"hydrogc"``).
        plip_complex_threshold : float
            Max distance (Å) for receptor atoms to include in interaction analysis.
        neighboring_residue_threshold : float
            Max distance (Å) for neighboring receptor residue detection.
        neighboring_ligand_threshold : float
            Max distance (Å) for neighboring ligand detection.
        data_dir : Path, optional
            Plinder data root for loading cofactors, affinity, etc.
        chain_to_seqres : dict[str, str], optional
            SEQRES per chain for binding affinity validation.
        ligand_smiles_dict : dict[str, str], optional
            Per-residue SMILES for components not in CCD (typically
            custom residues like Boltz's ``LIG``). When a residue's
            name appears in this dict, the user's SMILES takes
            precedence over CCD/PRD for both the canonical ``smiles``
            field and the stereo template used by
            :func:`_check_stereo_vs_template` — the caller is assumed
            to know that the CCD entry is absent or a placeholder.

        Returns
        -------
        Ligand or None
            Populated Ligand object, or None if no atoms found.
        """
        if data_dir is not None:
            global \
                COFACTORS, \
                ARTIFACTS, \
                LIST_OF_CCD_SYNONYMS, \
                CCD_SYNONYMS_DICT, \
                KINASE_INHIBITORS, \
                BINDING_AFFINITY
            if LIST_OF_CCD_SYNONYMS is None or CCD_SYNONYMS_DICT is None:
                LIST_OF_CCD_SYNONYMS, CCD_SYNONYMS_DICT = get_ccd_synonyms(data_dir)
            if COFACTORS is None:
                COFACTORS = parse_cofactors(data_dir)
            if ARTIFACTS is None:
                ARTIFACTS = parse_artifacts()
            if KINASE_INHIBITORS is None:
                KINASE_INHIBITORS = parse_kinase_inhibitors(data_dir)
            if BINDING_AFFINITY is None:
                try:
                    BINDING_AFFINITY = get_binding_affinity(data_dir)
                except Exception as e:
                    LOG.warning(f"Failed to load binding affinity data: {e}")
                    BINDING_AFFINITY = {"pchembl": {}, "target_sequence": {}}

        ligand_instance_chain = f"{ligand_instance}.{ligand_chain.asym_id}"

        # Select ligand atoms from biounit (AtomArray)
        lig_mask = (biounit.chain_id == ligand_instance_chain) & np.isin(
            biounit.res_id, residue_numbers
        )
        if not np.any(lig_mask):
            LOG.warning(f"from_pli: no ligand atoms for {ligand_instance_chain}")
            return None

        # Find complete residues within threshold distance of ligand
        lig_coords = biounit.coord[lig_mask]
        cell_list = struc.CellList(biounit, plip_complex_threshold)
        nearby_atom_mask = np.zeros(len(biounit), dtype=bool)
        for coord in lig_coords:
            indices = cell_list.get_atoms(coord, radius=plip_complex_threshold)
            nearby_atom_mask[indices[indices >= 0]] = True
        # Expand to complete residues to avoid broken aromatic rings
        nearby_mask = np.any(
            struc.get_residue_masks(biounit, np.where(nearby_atom_mask)[0]),
            axis=0,
        )
        nearby_atoms = biounit[nearby_mask]

        # Bonds propagate from biounit through array slicing;
        # only re-derive if missing
        if nearby_atoms.bonds is None:
            nearby_atoms.bonds = struc.connect_via_residue_names(nearby_atoms)

        # Split into receptor/ligand/water/metal
        receptor_mask = struc.filter_amino_acids(
            nearby_atoms
        ) | struc.filter_nucleotides(nearby_atoms)
        ligand_mask_local = nearby_atoms.chain_id == ligand_instance_chain
        water_mask = struc.filter_solvent(nearby_atoms)
        metal_mask = struc.filter_monoatomic_ions(nearby_atoms) & ~ligand_mask_local

        receptor_arr = nearby_atoms[receptor_mask & ~water_mask & ~metal_mask]
        ligand_arr = nearby_atoms[ligand_mask_local & ~water_mask]
        water_arr = nearby_atoms[water_mask]
        metal_arr = nearby_atoms[metal_mask]

        if receptor_arr.array_length() == 0 or ligand_arr.array_length() == 0:
            LOG.warning(
                f"from_pli: empty receptor or ligand for {ligand_instance_chain}"
            )
            return None

        # Chain mapping: chain_id is already in instance.asym format
        inv_mapping = {c: c for c in np.unique(nearby_atoms.chain_id)}

        peppr_interactions, peppr_waters = run_peppr_interactions(
            receptor_arr,
            ligand_arr,
            water_arr,
            metal_arr,
            ligand_instance_chain,
            inv_mapping,
        )

        # Get CCD codes from ligand atoms (one per residue, preserving duplicates)
        lig_atoms = biounit[lig_mask]
        ccd_code = "-".join(
            lig_atoms.res_name[lig_atoms.res_id == rn][0]
            for rn in residue_numbers
            if np.any(lig_atoms.res_id == rn)
        )
        # Get SMILES from CCD template via biotite, fall back to structure
        from plinder.core.structure.atoms import is_hydrogen_isotope

        smiles = None
        lig_heavy = lig_atoms[~is_hydrogen_isotope(lig_atoms.element)]
        res_names = list(
            dict.fromkeys(
                lig_heavy.res_name[lig_heavy.res_id == rn][0]
                for rn in residue_numbers
                if np.any(lig_heavy.res_id == rn)
            )
        )
        if len(res_names) == 1:
            resname = res_names[0]
            # User-supplied SMILES takes precedence — when the caller
            # explicitly provided one, CCD is assumed to be wrong or a
            # generic placeholder (biotite returns one for some codes
            # like "LIG"). Fall through to CCD then PRD otherwise.
            if ligand_smiles_dict and resname in ligand_smiles_dict:
                smiles = ligand_smiles_dict[resname]
            else:
                ccd_smiles = _get_ccd_smiles(resname)
                if ccd_smiles is None and resname.startswith("PRD_"):
                    ccd_smiles = _get_prd_smiles(resname)
                if ccd_smiles is not None:
                    smiles = ccd_smiles
        # Build per-residue custom stereo templates from user SMILES (only
        # populated for custom CIFs via from_custom_cif_file). The CIF atom
        # names for each residue are taken in file order, matching the
        # SMILES-parse-order assumption used for bond assignment.
        custom_templates: dict[str, Chem.Mol] | None = None
        if ligand_smiles_dict:
            custom_templates = {}
            for resname, user_smiles in ligand_smiles_dict.items():
                res_mask = lig_heavy.res_name == resname
                if not np.any(res_mask):
                    continue
                atom_names = list(lig_heavy.atom_name[res_mask])
                tmpl = _template_from_user_smiles(resname, user_smiles, atom_names)
                if tmpl is not None:
                    custom_templates[resname] = tmpl

        # Build the resolved (from 3D) mol once. It drives:
        #   - resolved_smiles (bond orders from CCD, stereo from 3D coords)
        #   - stereo match check against the CCD template (or custom SMILES)
        #   - fallback SMILES when the CCD/PRD/user-SMILES lookup failed
        resolved_smiles: str | None = None
        stereo_matches: bool | None = None
        try:
            from plinder.data.utils.annotations.cif_utils import atoms_to_rdkit_mol

            # biotite has no chiral tags → stereo assigned from 3D inside helper
            resolved_mol = atoms_to_rdkit_mol(lig_heavy)
            resolved_smiles = str(Chem.MolToSmiles(resolved_mol))
            # Compare resolved 3D stereo with CCD template stereo
            # (works for both single- and multi-residue ligands)
            stereo_matches = _check_stereo_vs_template(
                resolved_mol, custom_templates=custom_templates
            )
        except Exception as e:
            LOG.warning(f"Failed to compute resolved SMILES for {ccd_code}: {e}")
        # Fall back to resolved SMILES if no upstream source yielded one
        if smiles is None:
            smiles = resolved_smiles
        # Centroid
        centroid = list(lig_atoms.coord.mean(axis=0))
        ligand = cls(
            pdb_id=pdb_id,
            biounit_id=biounit_id,
            asym_id=ligand_chain.asym_id,
            instance=ligand_instance,
            ccd_code=ccd_code,
            plip_type=get_chain_type(ligand_chain.chain_type_str),
            bird_id=list(ligand_chain.mappings.get("BIRD", {"": None}))[0],  # type: ignore
            centroid=centroid,
            smiles=smiles or "",
            neighboring_residue_threshold=neighboring_residue_threshold,
            neighboring_ligand_threshold=neighboring_ligand_threshold,
            resolved_smiles=resolved_smiles or "",
            resolved_stereo_matches_template=stereo_matches,
            residue_numbers=residue_numbers,
        )

        # Find neighboring polymer residues (protein + nucleic acid) within threshold
        polymer_mask = struc.filter_amino_acids(biounit) | struc.filter_nucleotides(
            biounit
        )
        polymer_atoms = biounit[polymer_mask]
        if polymer_atoms.array_length() > 0:
            neighbor_cell = struc.CellList(
                polymer_atoms, ligand.neighboring_residue_threshold
            )
            near_poly_mask = np.zeros(len(polymer_atoms), dtype=bool)
            for coord in lig_coords:
                indices = neighbor_cell.get_atoms(
                    coord, radius=ligand.neighboring_residue_threshold
                )
                near_poly_mask[indices[indices >= 0]] = True
            near_prot = polymer_atoms[near_poly_mask]
        else:
            near_prot = polymer_atoms[:0]  # empty

        (
            ligand.num_neighboring_ppi_atoms_within_4A_of_gap,
            ligand.num_neighboring_ppi_atoms_within_8A_of_gap,
            ligand.num_missing_ppi_interface_residues,
            ligand.num_pli_atoms_within_4A_of_gap,
            ligand.num_pli_atoms_within_8A_of_gap,
            ligand.num_missing_pli_interface_residues,
        ) = annotate_interface_gaps_per_chain(
            interface_proximal_gaps, ligand_chain.asym_id
        )

        for chain_id in np.unique(near_prot.chain_id):
            if chain_id == ligand.instance_chain:
                continue
            # Skip chains classified as ligands — they belong in
            # neighboring_ligands/interacting_ligands, not neighboring_residues
            asym = chain_id.split(".")[-1] if "." in chain_id else chain_id
            if asym in ligand_like_chains:
                continue
            chain_atoms = near_prot[near_prot.chain_id == chain_id]
            resnums = list(dict.fromkeys(int(r) for r in chain_atoms.res_id))
            ligand.neighboring_residues[chain_id] = resnums
            # Store SEQRES for binding affinity validation
            asym_id = chain_id.split(".")[-1] if "." in chain_id else chain_id
            if chain_to_seqres and asym_id in chain_to_seqres:
                ligand.receptor_seqres[chain_id] = chain_to_seqres[asym_id]

        neighboring_asym_ids = {
            c.split(".")[-1]
            for c in np.unique(near_prot.chain_id)
            if c != ligand.instance_chain
        }

        ligand.covalent_linkages = extract_ligand_links_to_neighbouring_chains(
            all_covalent_dict, ligand.asym_id, neighboring_asym_ids, link_type="covale"
        )
        ligand.is_covalent = len(ligand.covalent_linkages) > 0

        # Find neighboring ligand chains
        near_lig_cell = struc.CellList(biounit, ligand.neighboring_ligand_threshold)
        near_lig_mask = np.zeros(len(biounit), dtype=bool)
        for coord in lig_coords:
            indices = near_lig_cell.get_atoms(
                coord, radius=ligand.neighboring_ligand_threshold
            )
            near_lig_mask[indices[indices >= 0]] = True
        near_all = biounit[near_lig_mask]

        ligand.neighboring_ligands = list(
            set(
                c
                for c in np.unique(near_all.chain_id)
                if c != ligand.instance_chain
                and "." in c
                and c.split(".")[1] in ligand_like_chains
            )
        )
        water_chains = set(
            c
            for c in np.unique(biounit.chain_id)
            if struc.filter_solvent(biounit[biounit.chain_id == c]).all()
        )
        # Populate interactions and waters from peppr results
        ligand.interactions = peppr_interactions
        ligand.waters = defaultdict(list)
        for w_chain, w_resnum in peppr_waters:
            ligand.waters[w_chain].append(w_resnum)

        # Derive interacting residues from peppr interaction hashes
        for instance_chain, residues in peppr_interactions.items():
            if instance_chain == ligand.instance_chain:
                continue
            if instance_chain in water_chains:
                continue
            if instance_chain.split(".")[1] in ligand_like_chains:
                ligand.interacting_ligands.append(instance_chain)
            else:
                if instance_chain not in ligand.interacting_residues:
                    ligand.interacting_residues[instance_chain] = []
                ligand.interacting_residues[instance_chain].extend(
                    int(r) for r in residues.keys()
                )
        # add rdkit properties and type assignments
        ligand.set_rdkit()
        if data_dir is not None:
            # set is_artifact and is_cofactor and is_other
            ligand.identify_artifacts_cofactors_and_other()
            # unique code parsing!
            ligand.unique_ccd_code = get_unique_ccd_longname(ligand.ccd_code)

        return ligand

    @cached_property
    def selection(self) -> str:
        """
        __Selection string for ligand
        """
        residue_selection = " or ".join(f"rnum={rnum}" for rnum in self.residue_numbers)
        ligand_selection = f"cname='{self.instance_chain}'"
        if len(self.residue_numbers):
            ligand_selection += f"and ({residue_selection})"
        return ligand_selection

    @cached_property
    def protein_chains_asym_id(self) -> list[str]:
        """Receptor chain IDs (protein/NA) within neighboring threshold of ligand.

        Returns empty list if the ligand is an artifact.
        """
        if self.is_artifact:
            return []
        else:
            return list(sorted(self.neighboring_residues.keys()))

    @cached_property
    def num_interacting_residues(self) -> int:
        """
        Number of residues interacting with a given ligand.
        """
        return sum(
            len(self.interacting_residues[chain]) for chain in self.interacting_residues
        )

    @cached_property
    def num_neighboring_residues(self) -> int:
        """Total count of receptor residues (protein/NA) within neighboring threshold."""
        return sum(
            len(self.neighboring_residues[chain]) for chain in self.neighboring_residues
        )

    @cached_property
    def is_proper(self) -> bool:
        """
        Check if ligand is a proper ligand (not an ion or artifact)
        """
        return not self.is_ion and not self.is_artifact

    @cached_property
    def num_interactions(self) -> int:
        """
        Number of interactions for a given ligand.
        """
        return sum(
            sum(len(i) for i in self.interactions[chain].values())
            for chain in self.interactions
        )

    @cached_property
    def num_unique_interactions(self) -> int:
        """
        Number of unique interactions
        """
        return sum(
            sum(len(set(i)) for i in self.interactions[chain].values())
            for chain in self.interactions
        )

    @cached_property
    def pocket_residues(self) -> dict[str, dict[int, str]]:
        """
        __Residues in the ligand's binding pocket which includes neighboring and interacting residues.
        """
        residues: dict[str, dict[int, str]] = {}
        for chain in self.neighboring_residues:
            if chain not in residues:
                residues[chain] = {}
            for residue in self.neighboring_residues[chain]:
                residues[chain][residue] = "neighboring"
        for chain in self.interacting_residues:
            if chain not in residues:
                residues[chain] = {}
            for residue in self.interacting_residues[chain]:
                residues[chain][residue] = "interacting"
        return residues

    def get_pocket_residues_set(self) -> dict[tuple[str, int], set[str]]:
        """
        Get a dict of pocket residues in the format (chain_id, residue_number)
        mapping to biounit instance set
        """
        pocket_residues_set = defaultdict(set)
        for chain in self.pocket_residues:
            for residue_number in self.pocket_residues[chain]:
                pocket_residues_set[(chain.split(".")[1], residue_number)].add(
                    chain.split(".")[0]
                )
        return pocket_residues_set

    def label_crystal_contacts(
        self,
        symmetry_mate_contacts: dict[
            tuple[str, int], dict[tuple[str, int], dict[int, set[int]]]
        ],
    ) -> None:
        """
        Label ligand contacts to chains that are not part of the biounit.
        """
        crystal_contacts: dict[tuple[str, int], set[int]] = defaultdict(set[int])

        # get contacts from neigchboring chain residues within the biounit
        pocket_residues = self.get_pocket_residues_set()

        for residue_number in self.residue_numbers:
            # get all inter-chain contacts for a given ligand
            contacts = symmetry_mate_contacts.get(
                (self.asym_id, residue_number), dict()
            )
            for x, y in contacts.items():
                # x is a tuple rec (chain_id, residue_number)
                # y is a dict of ligand atom_id : {image_idx} - set of symmetry operations
                num_crystal_image_contacts = len(y.values())
                # if detected contacts have more images than contact instances in the biounit pocket
                # then we assume that this is a crystal contact with a symmetry mate
                if num_crystal_image_contacts > len(pocket_residues.get(x, set())):
                    # on the edge cases it may not be clear which atom is in contact with the symmetry mate, thus better to store all?
                    for atom_id, image_idx in y.items():
                        crystal_contacts[x] |= {atom_id}
        # set crystal contacts
        self.crystal_contacts = crystal_contacts

    @cached_property
    def num_crystal_contacted_residues(self) -> int:
        """
        Number of residues from other symmetry mates which are in contact with this ligand.
        """
        return len(self.crystal_contacts)

    @cached_property
    def num_atoms_with_crystal_contacts(self) -> int:
        """
        Number of atoms in this ligand which are in contact with residues from other symmetry mates.
        """
        all_atoms = set()
        for x in self.crystal_contacts.values():
            all_atoms |= x
        return len(all_atoms)

    @cached_property
    def fraction_atoms_with_crystal_contacts(self) -> float | None:
        """
        Fraction of atoms in this ligand which are in contact with residues from other symmetry mates.
        """
        if self.num_heavy_atoms is None:
            return None
        return self.num_atoms_with_crystal_contacts / self.num_heavy_atoms

    @cached_property
    def num_pocket_residues(self) -> int:
        """
        Number of residues in the ligand's binding pocket.
        """
        return sum([len(self.pocket_residues[chain]) for chain in self.pocket_residues])

    @cached_property
    def id(self) -> str:
        """
        Unique identifier for a given ligand.
        """
        return "__".join([self.pdb_id, self.biounit_id, self.instance_chain])

    @cached_property
    def instance_chain(self) -> str:
        """
        Instance chain for a given ligand.
        """
        return f"{self.instance}.{self.asym_id}"

    @cached_property
    def interactions_counter(self) -> dict[str, dict[int, ty.Counter[str]]]:
        """
        __Counter of interactions for a given ligand.
        """
        interactions_counter: dict[str, dict[int, ty.Counter[str]]] = {}
        for chain in self.interactions:
            interactions_counter[chain] = {}
            for residue in self.interactions[chain]:
                interactions_counter[chain][residue] = Counter(
                    self.interactions[chain][residue]
                )
        return interactions_counter

    @cached_property
    def is_kinase_inhibitor(self) -> bool:
        """
        Check if ligand is a kinase inhibitor.
        """
        global KINASE_INHIBITORS
        if KINASE_INHIBITORS is None:
            data_dir = Path(get_config().data.plinder_dir)
            KINASE_INHIBITORS = parse_kinase_inhibitors(data_dir)
        return any(c in KINASE_INHIBITORS for c in self.ccd_code.split("-"))

    @cached_property
    def binding_affinity(self) -> float | None:
        """Binding affinity (pKd or pKi) from BindingDB when available.

        The affinity is only returned if the BindingDB target sequence
        matches at least one receptor chain SEQRES with 100% identity
        in the aligned core (terminal overhangs from tags/truncations
        are tolerated).  This guards against BindingDB's 85% sequence
        identity matching which can assign values to wrong complexes
        (see `#94 <https://github.com/plinder-org/plinder/issues/94>`_).
        """
        global BINDING_AFFINITY
        pdbid_ligid = f"{self.pdb_id}_{self.ccd_code}".upper()
        if BINDING_AFFINITY is None:
            data_dir = Path(get_config().data.plinder_dir)
            BINDING_AFFINITY = get_binding_affinity(data_dir)
        pchembl = BINDING_AFFINITY.get("pchembl", {})
        target_seqs = BINDING_AFFINITY.get("target_sequence", {})
        affinity = pchembl.get(pdbid_ligid)
        if affinity is None:
            return None
        # Validate: BindingDB target sequence must match a receptor chain
        bdb_seq = target_seqs.get(pdbid_ligid)
        if bdb_seq and self.receptor_seqres:
            if not any(
                sequences_match_core(bdb_seq, seq)
                for seq in self.receptor_seqres.values()
            ):
                LOG.warning(
                    f"binding_affinity: rejecting {pdbid_ligid} — "
                    "BindingDB target sequence does not match any receptor chain"
                )
                return None
        return float(affinity)

    def identify_artifacts_cofactors_and_other(self) -> None:
        """Set ``is_artifact``, ``is_cofactor``, and ``is_other`` flags in-place."""
        assert COFACTORS is not None
        assert ARTIFACTS is not None
        if self.ccd_code in COFACTORS:
            self.is_cofactor = True
        if self.ccd_code in ARTIFACTS:
            self.in_artifact_list = True

        if self.is_ion:
            self.is_artifact = False
        elif self.in_artifact_list:
            self.is_artifact = True
        elif lig_has_dummies(self.ccd_code):
            # check for dummy list including composites, too!
            self.is_artifact = True
        elif is_excluded_mol(self.smiles):
            self.is_artifact = True
        else:
            self.is_artifact = False

        # reset: artifacts should not count to other ligand class definitions
        if self.is_artifact:
            self.is_ion = False
            self.is_oligo = False
            self.is_cofactor = False
            self.is_lipinski = False
            self.is_fragment = False
            self.is_covalent = False

        # Indicator of whether a ligand type is not any of small molecule (lipinski, frag, coval), ion, cofactor, oligopeptide, oligosaccharide or oligopeptide.
        if not any(
            [
                self.is_invalid,
                self.is_ion,
                self.is_oligo,
                self.is_artifact,
                self.is_cofactor,
                self.is_lipinski,
                self.is_fragment,
                self.is_covalent,
            ]
        ):
            self.is_other = True

    def format_chains(
        self,
        chain_type: str,
        chains: dict[str, Chain],
    ) -> dict[str, ty.Any]:
        """
        Format chains for pd.DataFrame

        Parameters
        ----------
        self : Ligand
            Ligand object
        chain_type: str
            Chain tyoe
        chains: dict[str, Chain]
            Chain id : chain mapping
        Returns
        -------
        dict[str, str]
        """
        if chain_type == "protein":
            sub_chains = self.protein_chains_asym_id
        elif chain_type == "interacting_ligand":
            sub_chains = self.interacting_ligands
        elif chain_type == "neighboring_ligand":
            sub_chains = self.neighboring_ligands
        else:
            raise ValueError(f"chain_type={chain_type} not understood")
        sub_chains_data = [
            chains[instance_chain.split(".")[-1]].format(
                int(instance_chain.split(".")[0])
            )
            for instance_chain in sub_chains
        ]
        data: dict[str, list[ty.Any]] = defaultdict(list)
        if len(sub_chains_data) == 0:
            return {}
        for sub_chain in sub_chains_data:
            for key in sub_chain:
                data[f"ligand_{chain_type}_chains_{key}"].append(sub_chain[key])
        return data

    def format_residues(
        self, residue_type: str, chains: dict[str, Chain]
    ) -> dict[str, list[str]]:
        """
        Format residues for pd.DataFrame

        Parameters
        ----------
        self : Ligand
            Ligand object
        residue_type : str
            Chain tyoe
        chains : dict[str, Chain]
            Chain id : chain mapping

        Returns
        -------
        List of residues in the format "<chain>_<residue_number>_<residue_index>_<auth_number>"
        dict[str, list[str]]
        """
        if residue_type == "interacting":
            residues = self.interacting_residues
        elif residue_type == "neighboring":
            residues = self.neighboring_residues
        res = []
        for instance_chain in residues:
            _, chain = instance_chain.split(".")
            for residue_number in residues[instance_chain]:
                res.append(
                    f"{instance_chain}_{residue_number}_{chains[chain].residues[residue_number].index}_{chains[chain].residues[residue_number].auth_number}"
                )  # TODO: move some of this logic to Residue
        return {f"ligand_{residue_type}_residues": res}

    def format_interactions(self) -> dict[str, list[str]]:
        """
        Format interactions for pd.DataFrame

        Parameters
        ----------
        self : Ligand

        Returns
        -------
        List of interactions in the format "<chain>_<residue_number>_<interaction>"
        dict[str, list[str]]

        """
        interactions: list[str] = []
        for chain in self.interactions:
            for residue in self.interactions[chain]:
                for interaction in self.interactions[chain][int(residue)]:
                    interactions.append(f"{chain}_{residue}_{interaction}")
        return {"ligand_interactions": interactions}

    def format(self, chains: dict[str, Chain]) -> dict[str, ty.Any]:
        """Serialize ligand annotations to a flat dict for DataFrame export."""
        data: dict[str, ty.Any] = defaultdict(str)
        ignore_fields = set(
            [
                "posebusters_result",
                "interactions",
                "protein_chains",
                "interacting_ligands",
                "neighboring_ligands",
                "interacting_residues",
                "neighboring_residues",
                "pocket_residues",
            ]
        )
        for field, desc_type in self.get_descriptions_and_types().items():
            # blacklist fields that will be added with custom formatters below or that we don't want to add to the plindex
            descr = str(desc_type[0]).lstrip().replace("\n", " ")
            if descr.startswith("__") or field in ignore_fields:
                continue
            name = f"ligand_{field}"
            data[name] = getattr(self, field, None)

        # posebusters
        if self.posebusters_result is not None:
            for k in self.posebusters_result:
                data[f"ligand_posebusters_{k}"] = self.posebusters_result[k]
        # interactions
        data.update(self.format_interactions())
        # chains
        data.update(
            {"ligand_auth_id": chains[self.asym_id].auth_id}
        )  # not a cached_property b/c it needs chains!
        for chain_type in [
            "protein",
            "interacting_ligand",
            "neighboring_ligand",
        ]:
            data.update(self.format_chains(chain_type, chains))
        # residues
        for residue_type in ["interacting", "neighboring"]:
            data.update(self.format_residues(residue_type, chains))

        return data
