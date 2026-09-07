# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
from __future__ import annotations

import logging
import typing as ty
from collections import Counter, defaultdict
from dataclasses import dataclass
from functools import cache, cached_property
from pathlib import Path

import biotite.structure as struc
import numpy as np
import numpy.typing as npt
from pydantic import BeforeValidator, Field
from rdkit import Chem, RDLogger
from rdkit.Chem import QED, Crippen, rdMolDescriptors
from rdkit.Chem.rdchem import Mol

from plinder.core.utils.config import get_config
from plinder.core.utils.constants import BASE_DIR
from plinder.core.utils.sanitize import mol_from_smiles
from plinder.data.annotations.interaction_utils import (
    extract_ligand_links_to_neighbouring_chains,
    run_peppr_interactions,
)
from plinder.data.annotations.protein_utils import Chain, sequences_match_core
from plinder.data.annotations.utils import (
    DocBaseModel,
    description_excluded_from_flat_export,
)

LOG = logging.getLogger(__name__)


@dataclass
class BiounitSpatialIndex:
    """Reusable spatial and hierarchy index for one biological assembly."""

    cell_list: struc.CellList
    residue_starts: npt.NDArray[np.int_]
    chain_segments: dict[str, list[tuple[int, int]]]
    max_radius: float
    bond_neighbors: npt.NDArray[np.integer[ty.Any]] | None
    bond_types: npt.NDArray[np.integer[ty.Any]] | None

    @classmethod
    def from_atoms(
        cls,
        atoms: struc.AtomArray,
        max_radius: float,
    ) -> "BiounitSpatialIndex":
        """Build the expensive whole-assembly indexes exactly once."""
        if max_radius <= 0:
            raise ValueError("max_radius must be positive")
        chain_starts = struc.get_chain_starts(atoms, add_exclusive_stop=True)
        chain_segments: dict[str, list[tuple[int, int]]] = defaultdict(list)
        for start, stop in zip(chain_starts[:-1], chain_starts[1:]):
            chain_segments[str(atoms.chain_id[start])].append((int(start), int(stop)))
        if atoms.bonds is None:
            bond_neighbors = None
            bond_types = None
        else:
            bond_neighbors, bond_types = atoms.bonds.get_all_bonds()
        return cls(
            cell_list=struc.CellList(atoms, max_radius),
            residue_starts=ty.cast(
                npt.NDArray[np.int_],
                struc.get_residue_starts(atoms, add_exclusive_stop=True),
            ),
            chain_segments=dict(chain_segments),
            max_radius=max_radius,
            bond_neighbors=bond_neighbors,
            bond_types=bond_types,
        )

    @property
    def chain_ids(self) -> list[str]:
        """Return chain IDs in their original assembly order."""
        return list(self.chain_segments)

    def atom_indices_for_chain(self, chain_id: str) -> npt.NDArray[np.int_]:
        """Return sorted atom indices for a possibly segmented chain."""
        segments = self.chain_segments.get(chain_id, [])
        if not segments:
            return np.array([], dtype=int)
        return ty.cast(
            npt.NDArray[np.int_],
            np.concatenate(
                [np.arange(start, stop, dtype=int) for start, stop in segments]
            ),
        )

    def atom_indices_near(
        self,
        coordinates: npt.NDArray[np.floating[ty.Any]],
        radius: float,
    ) -> npt.NDArray[np.int_]:
        """Return sorted unique atom indices within *radius* of coordinates."""
        if radius > self.max_radius:
            raise ValueError(
                f"query radius {radius} exceeds indexed radius {self.max_radius}"
            )
        nearby: list[npt.NDArray[np.int_]] = []
        for coordinate in coordinates:
            indices = ty.cast(
                npt.NDArray[np.int_],
                self.cell_list.get_atoms(coordinate, radius=radius),
            )
            valid = indices[indices >= 0]
            if valid.size:
                nearby.append(valid)
        if not nearby:
            return np.array([], dtype=int)
        return ty.cast(npt.NDArray[np.int_], np.unique(np.concatenate(nearby)))

    def complete_residue_indices_near(
        self,
        coordinates: npt.NDArray[np.floating[ty.Any]],
        radius: float,
    ) -> npt.NDArray[np.int_]:
        """Expand nearby atom hits to their complete residues."""
        nearby = self.atom_indices_near(coordinates, radius)
        if nearby.size == 0:
            return nearby
        residue_indices = np.unique(
            np.searchsorted(self.residue_starts, nearby, side="right") - 1
        )
        return ty.cast(
            npt.NDArray[np.int_],
            np.concatenate(
                [
                    np.arange(
                        self.residue_starts[index],
                        self.residue_starts[index + 1],
                        dtype=int,
                    )
                    for index in residue_indices
                ]
            ),
        )

    def take_atoms(
        self,
        atoms: struc.AtomArray,
        indices: npt.NDArray[np.int_],
        *,
        include_bonds: bool,
    ) -> struc.AtomArray:
        """Take a subarray without globally reindexing every bond.

        ``AtomArray.__getitem__`` delegates to ``BondList.__getitem__``, which
        scans the full assembly bond graph even when only a few local atoms are
        requested. Rebuilding the selected adjacency from ``get_bonds()`` is
        equivalent and scales with the local slice instead.
        """
        subset = struc.AtomArray(len(indices))
        subset.coord = atoms.coord[indices]
        for category in atoms.get_annotation_categories():
            subset.set_annotation(category, atoms.get_annotation(category)[indices])
        if atoms.box is not None:
            subset.box = atoms.box.copy()
        if (
            include_bonds
            and self.bond_neighbors is not None
            and self.bond_types is not None
        ):
            selected = {
                int(global_index): local_index
                for local_index, global_index in enumerate(indices)
            }
            bonds = struc.BondList(len(indices))
            for local_index, global_index in enumerate(indices):
                neighbors = self.bond_neighbors[global_index]
                bond_types = self.bond_types[global_index]
                for neighbor, bond_type in zip(neighbors, bond_types):
                    if neighbor < 0:
                        continue
                    other_local_index = selected.get(int(neighbor))
                    if (
                        other_local_index is not None
                        and local_index < other_local_index
                    ):
                        bonds.add_bond(
                            local_index,
                            other_local_index,
                            int(bond_type),
                        )
            subset.bonds = bonds
        return subset


def get_water_chain_ids(atoms: struc.AtomArray) -> set[str]:
    """Return chain IDs whose atoms are all classified as solvent.

    Computing the solvent mask once avoids repeatedly slicing a potentially
    large biological assembly for every ligand in that assembly.
    """
    solvent_mask = struc.filter_solvent(atoms)
    solvent_chains = set(str(chain_id) for chain_id in atoms.chain_id[solvent_mask])
    non_solvent_chains = set(
        str(chain_id) for chain_id in atoms.chain_id[~solvent_mask]
    )
    return solvent_chains - non_solvent_chains


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
    # Tolerant parse (over-valent boron/main-group centres), consistent with the
    # rest of the ligand-SMILES pipeline; None (unparseable) falls back to no
    # stereo check, as before. mol_from_smiles perceives stereo (atom @/@@ and
    # double-bond E/Z), which compare_stereo_to_template reads.
    mol = mol_from_smiles(smiles)
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
    # Stereo comes from the SMILES parity tags (@/@@ and /\), which the parse
    # records on the atoms/bonds and RemoveHs(sanitize=False) preserves;
    # compare_stereo_to_template reads them to check chiral handedness and E/Z.
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
    ccd_code_dict: dict[str, str] | None = None,
) -> bool | None:
    """Compare resolved 3D stereo against a stereo template per residue.

    Template source precedence:
      1. ``custom_templates[resname]`` if provided — user-supplied SMILES
         templates win over CCD because the caller explicitly knows CCD
         is wrong or missing (biotite ships a generic placeholder for
         some codes like ``LIG`` that would otherwise silently hide
         stereo mismatches).
      2. A mapped CCD template from ``ccd_code_dict``, with its CCD atoms
         graph-matched and renamed to the custom CIF atom names.
      3. :func:`_get_ccd_mol(resname)` — CCD ideal coordinates.
      4. Return ``None`` for this residue if neither source yields a
         template.

    Delegates to :func:`compare_stereo_to_template` for the actual
    handedness comparison. Handles multi-residue ligands (e.g. glycans) by
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
        fragment_mol = frag.GetMol()

        # User-supplied custom templates take precedence over CCD: if
        # the caller provided a SMILES template for this residue, they
        # explicitly know CCD is wrong or missing (biotite ships a
        # generic placeholder for some codes like "LIG" that would
        # otherwise hide stereo mismatches).
        template_mol = None
        if custom_templates is not None:
            template_mol = custom_templates.get(resname)
        if template_mol is None:
            reference_code = (ccd_code_dict or {}).get(resname, resname)
            template_mol = _get_ccd_mol(reference_code)
            if template_mol is not None and reference_code != resname:
                template_mol = _template_with_fragment_atom_names(
                    template_mol,
                    fragment_mol,
                    custom_comp_id=resname,
                    reference_code=reference_code,
                )
        if template_mol is None:
            results.append(None)
            continue

        try:
            results.append(compare_stereo_to_template(fragment_mol, template_mol))
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


def _template_with_fragment_atom_names(
    template_mol: "Chem.Mol",
    fragment_mol: "Chem.Mol",
    *,
    custom_comp_id: str,
    reference_code: str,
) -> "Chem.Mol | None":
    """Rename a mapped CCD template using its graph match to a CIF residue."""
    if template_mol.GetNumAtoms() != fragment_mol.GetNumAtoms():
        LOG.warning(
            "Mapped CCD stereo template atom count differs for %s -> %s",
            custom_comp_id,
            reference_code,
        )
        return None
    match = fragment_mol.GetSubstructMatch(template_mol, useChirality=False)
    if len(match) != template_mol.GetNumAtoms():
        LOG.warning(
            "Mapped CCD stereo template cannot be graph-matched for %s -> %s",
            custom_comp_id,
            reference_code,
        )
        return None
    mapped = Chem.Mol(template_mol)
    for template_atom, fragment_index in zip(mapped.GetAtoms(), match):
        fragment_info = fragment_mol.GetAtomWithIdx(fragment_index).GetPDBResidueInfo()
        if fragment_info is None:
            return None
        info = Chem.AtomPDBResidueInfo()
        info.SetName(fragment_info.GetName())
        info.SetResidueName(custom_comp_id)
        info.SetResidueNumber(1)
        template_atom.SetMonomerInfo(info)
    return mapped


@cache
def _get_ccd_mol(comp_id: str) -> "Chem.Mol | None":
    """Return a sanitized RDKit Mol for a CCD component, or None if absent.

    Thin, None-safe wrapper over ``cif_utils._get_ccd_atomarray``. biotite has no
    stereo tags, so ``atoms_to_rdkit_mol`` perceives stereo (R/S and E/Z) from the
    ideal 3D coordinates via ``AssignStereochemistryFrom3D``.
    """
    from plinder.data.annotations.cif_utils import (
        _get_ccd_atomarray,
        atoms_to_rdkit_mol,
    )

    atoms = _get_ccd_atomarray(comp_id)
    if atoms is None:
        return None
    try:
        return atoms_to_rdkit_mol(atoms)
    except Exception as e:
        LOG.warning(f"Failed to build RDKit mol for {comp_id}: {e}")
        return None


def _get_ccd_smiles(comp_id: str) -> str | None:
    """Get SMILES from a CCD component, with stereochemistry from ideal 3D."""
    mol = _get_ccd_mol(comp_id)
    if mol is None:
        return None
    return str(Chem.MolToSmiles(mol))


def lig_has_dummies(
    ligand_code: str,
    dummy_lig_list: list[str] = [
        # "DUM",  # OBS -> UNX (already listed below)
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
        Defaults to ['UNX', 'UNL', 'UNK', 'UPL', 'DN', 'N'].

    Returns:
        bool: if ligand considered as dummy and treated as artifact
    """
    # check for dummy list including composites, too!
    return len(set(ligand_code.split("-")).intersection(dummy_lig_list)) > 0


# lazy evaluate data fetches referenced as module globals
COFACTORS: set[str] | None = None
# RDKit canonical SMILES of the cofactor / artifact reference molecules. Matching
# on structure (not CCD code) classifies a ligand that is the same molecule under
# a *different current* code — the real "synonym" case. No "obsolete" codes.
COFACTOR_SMILES: set[str] | None = None
ARTIFACT_SMILES: set[str] | None = None
# instantiate artifact list once and reuse variable
ARTIFACTS: set[str] | None = None
BINDING_AFFINITY: dict[str, ty.Any] | None = None


_OLIGO_SMARTS = {
    "peptide": Chem.MolFromSmarts("C(=O)C[N;D2,D3]C(=[O;D1])CN"),
    "saccharide": Chem.MolFromSmarts("O-[C;R0,R1]-[C;R0,R1]-[O,S;R0;D2]-[C;R1]-[O;R1]"),
    "nucleotide": Chem.MolFromSmarts("P(=O)([O-,OH])(OC[C;r5])O[C;r5]"),
}
_POLYMER_CLASS_FIELDS = tuple(
    f"is_{size}{family}"
    for family in ("saccharide", "nucleotide", "peptide")
    for size in ("mono", "oligo")
)


def _has_nucleobase_attachment(mol: Mol, ring: set[int]) -> bool:
    """Return whether a sugar-ring carbon is bonded to a ring nitrogen."""
    for atom_index in ring:
        atom = mol.GetAtomWithIdx(atom_index)
        if atom.GetAtomicNum() != 6:
            continue
        for neighbor in atom.GetNeighbors():
            if (
                neighbor.GetIdx() not in ring
                and neighbor.GetAtomicNum() == 7
                and neighbor.IsInRing()
            ):
                return True
    return False


def _sugar_and_nucleotide_unit_counts(mol: Mol) -> tuple[int, int]:
    """Count conservative sugar and phosphorylated nucleoside ring motifs."""
    saccharide_units = 0
    nucleotide_units = 0
    phosphorus = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15]
    for atom_ring in mol.GetRingInfo().AtomRings():
        if len(atom_ring) not in {5, 6}:
            continue
        ring = set(atom_ring)
        atoms = [mol.GetAtomWithIdx(index) for index in atom_ring]
        atomic_numbers = [atom.GetAtomicNum() for atom in atoms]
        if any(atom.GetIsAromatic() for atom in atoms):
            continue
        if atomic_numbers.count(8) != 1 or atomic_numbers.count(6) != len(atoms) - 1:
            continue

        has_nucleobase = _has_nucleobase_attachment(mol, ring)
        close_to_phosphate = any(
            len(Chem.GetShortestPath(mol, ring_index, phosphorus_index)) - 1 <= 3
            for ring_index in ring
            for phosphorus_index in phosphorus
        )
        if has_nucleobase and close_to_phosphate:
            nucleotide_units += 1
            continue
        if has_nucleobase:
            # A nucleoside sugar is not a saccharide ligand class.
            continue

        external_oxygen_neighbors = {
            neighbor.GetIdx()
            for atom_index in ring
            for neighbor in mol.GetAtomWithIdx(atom_index).GetNeighbors()
            if neighbor.GetIdx() not in ring and neighbor.GetAtomicNum() == 8
        }
        if len(external_oxygen_neighbors) >= 2:
            saccharide_units += 1
    return saccharide_units, nucleotide_units


def _peptide_unit_count(mol: Mol) -> int:
    """Count alpha-amino-acid backbone centers in a molecule."""
    units = 0
    for atom in mol.GetAtoms():
        if (
            atom.GetAtomicNum() != 6
            or atom.GetHybridization() != Chem.HybridizationType.SP3
        ):
            continue
        neighbors = list(atom.GetNeighbors())
        if not any(neighbor.GetAtomicNum() == 7 for neighbor in neighbors):
            continue
        has_carbonyl = False
        for neighbor in neighbors:
            if neighbor.GetAtomicNum() != 6:
                continue
            for bond in neighbor.GetBonds():
                other = bond.GetOtherAtom(neighbor)
                if (
                    other.GetAtomicNum() == 8
                    and bond.GetBondType() == Chem.BondType.DOUBLE
                ):
                    has_carbonyl = True
                    break
            if has_carbonyl:
                break
        units += int(has_carbonyl)
    return units


@cache
def _smiles_descriptors(smiles: str | None) -> dict[str, ty.Any] | None:
    """Cached RDKit scalar descriptors for one canonical SMILES.

    Descriptors are a pure function of the SMILES, so they are memoized here
    (keyed on the string) instead of recomputed per ligand instance in
    ``set_rdkit``. Parsing goes through :func:`~plinder.core.utils.sanitize.mol_from_smiles` so
    over-valent ligands still yield descriptors rather than a silent ``None``.
    """
    mol = mol_from_smiles(smiles)
    if mol is None:
        return None
    descriptors: dict[str, ty.Any] = {
        "molecular_weight": rdMolDescriptors.CalcExactMolWt(mol),
        "num_rot_bonds": rdMolDescriptors.CalcNumRotatableBonds(mol),
        "num_hba": rdMolDescriptors.CalcNumHBA(mol),
        "num_hbd": rdMolDescriptors.CalcNumHBD(mol),
        "crippen_clogp": Crippen.MolLogP(mol),
        "num_rings": rdMolDescriptors.CalcNumRings(mol),
        "num_heavy_atoms": rdMolDescriptors.CalcNumHeavyAtoms(mol),
        "tpsa": rdMolDescriptors.CalcTPSA(mol),
    }
    # QED internally re-runs RemoveHs/sanitization, which rejects the over-valent
    # main-group centres (boron cages, hypervalent metals) that peppr.sanitize
    # tolerates. Keep the other descriptors and leave qed unset in that case.
    try:
        descriptors["qed"] = QED.qed(mol)
    except Exception:
        descriptors["qed"] = None
    return descriptors


@cache
def _classify_ligand_polymer_classes(
    smiles: str,
    sum_disconnected_units: bool = False,
) -> tuple[bool, ...]:
    """Return cached class flags for one SMILES representation.

    Disconnected fragments are normally treated as a mixture, so two copies
    of a monomer do not become an oligomer.  The caller may sum their units
    when separate residue records are known to belong to one multi-residue
    Parsing goes through :func:`~plinder.core.utils.sanitize.mol_from_smiles` so
    over-valent ligands are still classified rather than silently skipped.
    ligand whose inter-residue bonds were absent from the deposited CIF.
    """
    mol = mol_from_smiles(smiles) if smiles else None
    unit_counts = {"saccharide": 0, "nucleotide": 0, "peptide": 0}
    oligo_matches = {family: False for family in unit_counts}
    if mol is not None:
        for fragment in Chem.GetMolFrags(mol, asMols=True):
            saccharide_units, nucleotide_units = _sugar_and_nucleotide_unit_counts(
                fragment
            )
            peptide_units = _peptide_unit_count(fragment)
            if sum_disconnected_units:
                unit_counts["saccharide"] += saccharide_units
                unit_counts["nucleotide"] += nucleotide_units
                unit_counts["peptide"] += peptide_units
            else:
                unit_counts["saccharide"] = max(
                    unit_counts["saccharide"], saccharide_units
                )
                unit_counts["nucleotide"] = max(
                    unit_counts["nucleotide"], nucleotide_units
                )
                unit_counts["peptide"] = max(unit_counts["peptide"], peptide_units)
        for family, pattern in _OLIGO_SMARTS.items():
            if pattern is not None:
                oligo_matches[family] = mol.HasSubstructMatch(pattern)

    result: list[bool] = []
    for family in ("saccharide", "nucleotide", "peptide"):
        is_oligo = oligo_matches[family] or unit_counts[family] >= 2
        result.extend((unit_counts[family] == 1 and not is_oligo, is_oligo))
    return tuple(result)


def classify_ligand_polymer_classes(
    smiles: str | None,
    *,
    resolved_smiles: str | None = None,
    is_multi_residue: bool = False,
) -> dict[str, bool]:
    """Classify mono/oligo saccharide, nucleotide, and peptide ligands.

    Existing oligo SMARTS are retained. Conservative structural unit counts
    distinguish one unit from multiple units.  Resolved coordinates may add
    evidence that the canonical identity omitted; disconnected resolved units
    are summed only for a ligand already known to span multiple residues.
    """
    flags = dict(
        zip(
            _POLYMER_CLASS_FIELDS,
            _classify_ligand_polymer_classes(smiles or ""),
        )
    )
    if not resolved_smiles:
        return flags

    resolved_flags = dict(
        zip(
            _POLYMER_CLASS_FIELDS,
            _classify_ligand_polymer_classes(
                resolved_smiles,
                sum_disconnected_units=is_multi_residue,
            ),
        )
    )
    for family in ("saccharide", "nucleotide", "peptide"):
        mono_field = f"is_mono{family}"
        oligo_field = f"is_oligo{family}"
        is_oligo = flags[oligo_field] or resolved_flags[oligo_field]
        flags[oligo_field] = is_oligo
        flags[mono_field] = not is_oligo and (
            flags[mono_field] or resolved_flags[mono_field]
        )
    return flags


def _choose_ligand_smiles_by_heavy_atom_count(
    reference_smiles: str | None,
    resolved_smiles: str | None,
) -> str:
    """Choose the chemically useful identity, then the more complete one.

    A connected coordinate-derived molecule wins over disconnected
    per-component CCD templates because only it encodes observed
    inter-residue bonds. Otherwise, choose the representation containing more
    heavy atoms and prefer the resolved representation on a tie.
    """

    def valid_candidate(smiles: str | None) -> tuple[str, int, int] | None:
        if not smiles:
            return None
        mol = mol_from_smiles(smiles)
        if mol is None:
            return None
        return smiles, int(mol.GetNumHeavyAtoms()), len(Chem.GetMolFrags(mol))

    reference = valid_candidate(reference_smiles)
    resolved = valid_candidate(resolved_smiles)
    if reference is None:
        return resolved[0] if resolved is not None else (reference_smiles or "")
    if resolved is None:
        return reference[0]
    if resolved[2] == 1 and reference[2] > 1:
        return resolved[0]
    return reference[0] if reference[1] > resolved[1] else resolved[0]


def get_artifact_codes(data_dir: Path) -> set[str]:
    """Load the artifact CCD set needed for cheap ingest preflight."""
    global ARTIFACTS, ARTIFACT_SMILES

    artifacts = ARTIFACTS
    if artifacts is None:
        artifacts = parse_artifacts()
        ARTIFACTS = artifacts
        ARTIFACT_SMILES = _reference_smiles(artifacts)
    return artifacts


def _reference_smiles(codes: set[str]) -> set[str]:
    """RDKit canonical SMILES for a set of CCD codes (skipping unresolved ones).

    Classifying cofactors / artifacts by this SMILES set (in addition to the CCD
    code) catches a ligand that is the *same molecule* under a different current
    CCD code — the real "synonym" case. Uses :func:`_get_ccd_smiles` (derived from
    the CCD atoms), so it needs no external descriptor column, and the
    strings are directly comparable to a ligand's own ``self.smiles``.
    """
    smiles: set[str] = set()
    for code in codes:
        smi = _get_ccd_smiles(code)
        if smi:
            smiles.add(smi)
    return smiles


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
            # "CPG",  # OBS -> 5GP (already listed)
            # "G25",  # OBS, no REL successor
            "5GP",
        ],  # + ["GNP", "GTN"],  # GNP/GTN is inhibitor
        "Cytidine nucleotides": [
            "C",
            "C5P",
            "CDP",
            "CTP",
        ],  # C25 (OBS) -> C5P, already listed
        "Thymidine nucleotides": [
            "TMP",
            "DT",
            "TTP",
            "THM",
            "TYD",
        ],  # T (OBS) -> DT, already listed
        "Uridine nucleotides": [
            "U",
            "DU",
            "U5P",
            "UMP",
            "UDP",
            "UTP",
        ],  # U25 (OBS) -> U5P, already listed
        # "MIO": ["CRW"],  # CRW is OBS with no REL successor (kept; verify)
        "NAD": ["NAD"],  # was NAH (OBS -> NAD)
        "Glutathione": ["GPR"],  # was CYP (OBS -> GPR)
        "Biopterin": ["HBI", "H4B"],  # was HBL->HBI, BH4->H4B, THB->H4B (dup)
        "Tetrahydrofolic acid": ["MEF"],
        "Lumazine": ["DLZ"],
        "Menaquinone": ["MQ8", "MQ9", "MQE", "MQ7"],  # 7MQ (OBS) -> MQ7
        "Heme": ["1CP", "CP3", "MMP", "UP2", "UP3"],
        "Methanopterin": ["H4M", "H4Z"],
        "Lipoamide": ["LPM"],
        "Ubiquinone": ["DCQ", "HQE", "PLQ"],
        "Pyridoxal": ["PXL", "UEG"],
        "Siderophores": ["488", "EB4", "SE8"],
        "Methanofuran": ["MFN"],
        "Vitamin A": ["BCR", "ECH", "EQ3"],  # RAW (OBS) -> ECH, already listed
        "Vitamin K1": ["PQN"],
        "CHLOROPHYLL and similar": [
            "CLA",
            "CHL",
            "CL0",
            # "CL1",  # OBS, no REL successor
            # "CL2",  # OBS, no REL successor
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

    return cofactors


@cache
def parse_artifacts() -> set[str]:
    """Get and parse artifacts
    Returns:
        set[str]: set[str]
    """
    artifact_log = BASE_DIR / "annotations/static_files/artifacts_badlist.csv"
    with open(artifact_log, "r") as f:
        lines = f.readlines()
    artifacts = {l.strip() for l in lines if not l.startswith("#")}
    return artifacts


@cache
def get_binding_affinity(data_dir: Path) -> ty.Any:
    """Load BindingDB affinity data (pchembl values + target sequences)."""
    from plinder.data.pipeline.io import download_affinity_data

    return download_affinity_data(data_dir=data_dir)


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
    max_charge: int = 6,
    max_linear_hydrocarbon_linker: int = 12,
) -> bool:
    """Exclude some molecules by default as useless for druglikeness
    Uses OR logic for violating rules:
        - less than 2 carbon atoms
        - less than 5 non-hydrogen atoms
        - charge larger than +/- 6
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


def is_known_artifact_ligand(
    residue_names: ty.Iterable[str],
    artifact_codes: set[str],
) -> bool:
    """Return whether a ligand chain is provably an artifact before assembly.

    Composite ligands are only rejected from CCD-code rules that are identical
    to :meth:`Ligand.identify_artifacts_cofactors_and_other`.  The molecular
    exclusion rules are evaluated only for a single-residue CCD ligand, where
    the CCD template is the same canonical source used by ``Ligand.from_pli``.
    Uncertain cases deliberately return ``False`` and follow the full path.
    """
    names = [str(name) for name in residue_names]
    if not names:
        return False
    ccd_code = "-".join(names)
    if ccd_code in artifact_codes or lig_has_dummies(ccd_code):
        return True
    if len(names) != 1:
        return False
    smiles = _get_ccd_smiles(names[0])
    if smiles is None:
        return False
    try:
        return is_excluded_mol(smiles)
    except Exception as exc:
        LOG.warning("Could not preclassify CCD %s: %s", names[0], exc)
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
        description="[EXCLUDE] RCSB PDB ID, see https://mmcif.wwpdb.org/dictionaries/mmcif_pdbx_v50.dic/Items/_entry.id.html",
    )
    biounit_id: str = Field(
        default_factory=str,
        description="[EXCLUDE] Biounit id",
    )
    id_legacy: str = Field(
        default="",
        description="Historical ligand ID using global assembly-operation chain instances",
    )
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
    bird_id: str = Field(default_factory=str, description="Ligand BIRD (PRD) id")
    centroid: list[float] = Field(
        default_factory=list, description="Ligand center of geometry"
    )
    smiles: str = Field(
        default_factory=str,
        description="Ligand SMILES from CCD lookup (user-supplied SMILES wins for custom CIFs) or resolved 3D; for composite ligands, the valid representation with more heavy atoms is used and resolved connectivity wins ties",
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
        default_factory=list,
        description="[EXCLUDE] Ligand residue numbers",
    )
    member_residue_numbers: dict[str, list[int]] = Field(
        default_factory=dict,
        description="[EXCLUDE] Residue numbers per member instance-chain. A ligand may "
        "span several covalently-linked chains (e.g. a macrocycle whose parts "
        "are deposited as separate chains); this maps each member "
        "'{instance}.{asym_id}' to its residue numbers so every atom in the "
        "ligand can be selected. Single-chain ligands map their one instance-chain.",
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
        description="[CUSTOM_EXPORT] Dictionary of neighboring residues, with {instance}.{chain} key and residue number value",
    )
    neighboring_ligands: list[str] = Field(
        default_factory=list,
        description="[EXCLUDE] List of neighboring ligands {instance}.{chain}",
    )
    receptor_seqres: dict[str, str] = Field(
        default_factory=dict,
        description="[EXCLUDE] SEQRES sequences of neighboring receptor chains for affinity validation",
    )
    interacting_residues: dict[str, list[int]] = Field(
        default_factory=dict,
        description="[CUSTOM_EXPORT] Dictionary of interacting residues, with {instance}.{chain} key and residue number value",
    )
    interacting_ligands: list[str] = Field(
        default_factory=list,
        description="[EXCLUDE] List of interacting ligands {instance}.{chain}",
    )
    # TODO: rename interactions description; hash format kept for backward compatibility
    # (now computed by peppr, not PLIP)
    interactions: dict[str, dict[int, list[str]]] = Field(
        default_factory=dict,
        description="[EXCLUDE] Dictionary of {instance}.{chain} to residue number to list of interaction hashes",
    )

    @classmethod
    def document_properties(
        cls, prefix: str
    ) -> ty.Generator[tuple[str, str | None, str], ty.Any, ty.Any]:
        """Describe model fields plus the flat columns emitted by ``format()``."""
        yield from super().document_properties(prefix)
        custom_columns = (
            (
                "residue_numbers",
                "list[int]",
                "Resolved ligand residue numbers used to reconstruct this ligand "
                "from the source mmCIF",
            ),
            (
                "instance_chains",
                "list[str]",
                "All instance chains this ligand spans (several for a multi-chain "
                "covalent ligand); used to reconstruct the whole-molecule SDF from "
                "the source mmCIF",
            ),
            (
                "water_residues",
                "list[str]",
                "Interacting water residues encoded as "
                "<instance>.<asym>_<residue_number>",
            ),
            (
                "interactions",
                "list[str]",
                "Protein-ligand interactions encoded by receptor chain, residue "
                "number, and interaction type",
            ),
            ("auth_id", "str", "Author chain ID of the ligand"),
        )
        for suffix, dtype, description in custom_columns:
            yield f"{prefix}_{suffix}", dtype, description

    neighboring_residue_threshold: float = Field(
        default=6.0,
        description="[EXCLUDE] Maximum distance to consider receptor residues (protein/NA) neighboring",
    )
    neighboring_ligand_threshold: float = Field(
        default=4.0,
        description="[EXCLUDE] Maximum distance to consider ligands neighboring",
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
    is_monosaccharide: bool = Field(
        default=False,
        description="Indicator of whether a ligand contains one saccharide unit",
    )
    is_oligosaccharide: bool = Field(
        default=False,
        description="Indicator of whether a ligand contains multiple saccharide units",
    )
    is_mononucleotide: bool = Field(
        default=False,
        description="Indicator of whether a ligand contains one nucleotide unit",
    )
    is_oligonucleotide: bool = Field(
        default=False,
        description="Indicator of whether a ligand contains multiple nucleotide units",
    )
    is_monopeptide: bool = Field(
        default=False,
        description="Indicator of whether a ligand contains one peptide unit",
    )
    is_oligopeptide: bool = Field(
        default=False,
        description="Indicator of whether a ligand contains multiple peptide units",
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
        + "(Lipinski, Fragment or covalent), ion, cofactor, mono/oligo (peptide, saccharide or nucleotide) or artifact",
    )
    is_invalid: bool = Field(
        default=False, description="Indicator of whether a ligand is invalid"
    )
    unique_ccd_code: str | None = Field(
        default=None, description="Ligand representative CCD code after de-duplicating"
    )
    crystal_contacts: CrystalContacts = Field(
        default_factory=dict,
        description="[EXCLUDE] Dictionary of {chain} to residue number to set of interacting crystal contacts",
    )
    waters: dict[str, list[int]] = Field(
        default_factory=dict,
        description="[EXCLUDE] Dictionary of {instance}.{chain} to list of interacting water residue numbers",
    )
    """Ligand annotation dataclass.

    Holds structural, chemical, and interaction annotations for a single
    ligand chain in a protein–ligand (or NA–ligand) complex.
    """

    def set_rdkit(self) -> None:
        """Compute RDKit molecular descriptors from ``self.smiles``."""
        try:
            is_multi_residue = self._is_multi_residue
            if is_multi_residue:
                self.smiles = _choose_ligand_smiles_by_heavy_atom_count(
                    self.smiles,
                    self.resolved_smiles,
                )
            # Descriptors are a pure function of the SMILES, so they are looked up
            # from the SMILES-keyed cache rather than recomputed per instance.
            descriptors = _smiles_descriptors(self.smiles)
            if descriptors is not None:
                self.molecular_weight = descriptors["molecular_weight"]
                self.num_rot_bonds = descriptors["num_rot_bonds"]
                self.num_hba = descriptors["num_hba"]
                self.num_hbd = descriptors["num_hbd"]
                self.crippen_clogp = descriptors["crippen_clogp"]
                self.num_rings = descriptors["num_rings"]
                self.tpsa = descriptors["tpsa"]
                self.qed = descriptors["qed"]
                # The composite-ligand choice above may switch ``smiles`` to the
                # resolved representation, so keep the heavy-atom count (set at
                # construction from the reference SMILES) in step with it.
                self.num_heavy_atoms = descriptors["num_heavy_atoms"]
                if self.num_heavy_atoms and self.num_resolved_heavy_atoms:
                    self.num_unresolved_heavy_atoms = (
                        self.num_heavy_atoms - self.num_resolved_heavy_atoms
                    )
            # classify ligand based on the (tolerantly sanitized) molecule; skip
            # when the SMILES is empty/unparseable so a structurally valid ligand
            # that merely lost its SMILES (e.g. a multi-residue peptide) is not
            # flagged invalid.
            mol = mol_from_smiles(self.smiles)
            if mol is not None:
                self.classify_ligand_type(mol)

        except Exception as e:
            logging.warning(f"Error in setting rdkit for {self.id}: {e}")
            # Multi-residue ligands (peptides) may fail SMILES derivation
            # but are still structurally valid
            if self.smiles is None:
                self.is_invalid = True

    def classify_ligand_type(self, mol: Mol) -> None:
        """Classify ligand as ion, Lipinski, fragment, or mono/oligo class.

        Uses SMARTS patterns and Lipinski rules to assign granular type
        beyond the chain-type classification.

        Note
        ----
        Oligo smarts obtained from https://doi.org/10.1021/acs.jcim.3c01573

        Args:
            mol (Mol): RDKit compatible molecule
        """
        RDLogger.DisableLog("rdApp.*")
        if mol is not None:
            if is_single_atom_or_ion(mol):
                self.is_ion = True
            else:
                try:
                    polymer_classes = classify_ligand_polymer_classes(
                        self.smiles,
                        resolved_smiles=self.resolved_smiles,
                        is_multi_residue=self._is_multi_residue,
                    )
                    for field, value in polymer_classes.items():
                        setattr(self, field, value)
                    if self.plip_type == "SACCHARIDE" and self._is_multi_residue:
                        self.is_monosaccharide = False
                        self.is_oligosaccharide = True
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
        all_covalent_dict: dict[str, list[tuple[str, str]]],
        # TODO: rename plip_complex_threshold -> complex_threshold
        plip_complex_threshold: float = 10.0,
        neighboring_residue_threshold: float = 6.0,
        neighboring_ligand_threshold: float = 4.0,
        data_dir: ty.Optional[Path] = None,
        chain_to_seqres: dict[str, str] | None = None,
        ligand_smiles_dict: dict[str, str] | None = None,
        ligand_ccd_code_dict: dict[str, str] | None = None,
        water_chains: set[str] | None = None,
        spatial_index: BiounitSpatialIndex | None = None,
        member_residue_numbers: dict[str, list[int]] | None = None,
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
        all_covalent_dict : dict[str, list[tuple[str, str]]]
            Covalent linkages by type (``"covale"``, ``"metalc"``, ``"hydrog"``).
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
            precedence over CCD for both the canonical ``smiles``
            field and the stereo template used by
            :func:`_check_stereo_vs_template` — the caller is assumed
            to know that the CCD entry is absent or a placeholder.
        ligand_ccd_code_dict : dict[str, str], optional
            Custom component ID to reference CCD code. This changes the
            reported ``ccd_code`` while ``ligand_smiles_dict`` carries the
            resolved reference SMILES used for RDKit and stereo checks.
        water_chains : set[str], optional
            Chain IDs containing only solvent atoms. Pass a precomputed set
            when processing multiple ligands from the same assembly.
        spatial_index : BiounitSpatialIndex, optional
            Reusable whole-assembly spatial and hierarchy index.
        member_residue_numbers : dict[str, list[int]], optional
            Residue numbers per member instance-chain for a ligand that spans
            several covalently-linked chains (a macrocycle deposited as
            separate chains). Defaults to the single primary chain
            (``{ligand_instance_chain: residue_numbers}``) when omitted.

        Returns
        -------
        Ligand or None
            Populated Ligand object, or None if no atoms found.
        """
        if data_dir is not None:
            global COFACTORS, ARTIFACTS, COFACTOR_SMILES, ARTIFACT_SMILES
            global BINDING_AFFINITY
            if COFACTORS is None:
                COFACTORS = parse_cofactors(data_dir)
                COFACTOR_SMILES = _reference_smiles(COFACTORS)
            if ARTIFACTS is None:
                ARTIFACTS = parse_artifacts()
                ARTIFACT_SMILES = _reference_smiles(ARTIFACTS)
            if BINDING_AFFINITY is None:
                try:
                    BINDING_AFFINITY = get_binding_affinity(data_dir)
                except Exception as e:
                    LOG.warning(f"Failed to load binding affinity data: {e}")
                    BINDING_AFFINITY = {"pchembl": {}, "target_sequence": {}}

        ligand_instance_chain = f"{ligand_instance}.{ligand_chain.asym_id}"
        # A ligand may span several covalently-linked chains. Default to the
        # single primary chain when no members were supplied.
        if member_residue_numbers is None:
            member_residue_numbers = {ligand_instance_chain: residue_numbers}
        member_instance_chains = set(member_residue_numbers)
        member_asym_ids = {ic.split(".")[-1] for ic in member_instance_chains}

        if spatial_index is None:
            spatial_index = BiounitSpatialIndex.from_atoms(
                biounit,
                max(
                    plip_complex_threshold,
                    neighboring_residue_threshold,
                    neighboring_ligand_threshold,
                ),
            )

        # Select ligand atoms across all covalently-linked member chains
        # without scanning the full assembly.
        lig_index_parts = []
        for member_chain, member_rns in member_residue_numbers.items():
            member_indices = spatial_index.atom_indices_for_chain(member_chain)
            lig_index_parts.append(
                member_indices[np.isin(biounit.res_id[member_indices], member_rns)]
            )
        lig_indices = (
            np.concatenate(lig_index_parts)
            if lig_index_parts
            else np.array([], dtype=int)
        )
        if lig_indices.size == 0:
            LOG.warning(f"from_pli: no ligand atoms for {ligand_instance_chain}")
            return None
        lig_atoms = spatial_index.take_atoms(biounit, lig_indices, include_bonds=True)
        lig_coords = lig_atoms.coord

        # Find complete residues within threshold distance of ligand. Residue
        # starts and the whole-assembly CellList are shared by every ligand.
        nearby_indices = spatial_index.complete_residue_indices_near(
            lig_coords, plip_complex_threshold
        )
        nearby_atoms = spatial_index.take_atoms(
            biounit, nearby_indices, include_bonds=True
        )

        # build_biounit guarantees bonds and they propagate through the slice;
        # if not - fail loud!
        if nearby_atoms.bonds is None:
            raise ValueError(
                f"from_pli: pocket atoms for {ligand_instance_chain} arrived "
                "without bonds; build_biounit must supply them."
            )

        # Split into receptor/ligand/water/metal
        receptor_mask = struc.filter_amino_acids(
            nearby_atoms
        ) | struc.filter_nucleotides(nearby_atoms)
        ligand_mask_local = np.isin(nearby_atoms.chain_id, list(member_instance_chains))
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

        # CCD codes, one per residue in atom order across all member chains.
        # Keying on (chain, res_id) avoids collisions when merged chains reuse
        # residue numbers.
        def _residues_in_order(atoms: "struc.AtomArray") -> list[str]:
            seen: set[tuple[str, int]] = set()
            names: list[str] = []
            for chain_id, res_id, res_name in zip(
                atoms.chain_id, atoms.res_id, atoms.res_name
            ):
                key = (str(chain_id), int(res_id))
                if key not in seen:
                    seen.add(key)
                    names.append(str(res_name))
            return names

        residue_component_ids = _residues_in_order(lig_atoms)
        ccd_code = "-".join(
            (ligand_ccd_code_dict or {}).get(component_id, component_id)
            for component_id in residue_component_ids
        )
        # Get SMILES from CCD template via biotite, fall back to structure
        from biotite.structure import filter_heavy

        smiles = None
        lig_heavy = lig_atoms[filter_heavy(lig_atoms)]
        res_names = _residues_in_order(lig_heavy)
        reference_fragments: list[str] = []
        for resname in res_names:
            component_smiles: str | None
            # User-supplied SMILES takes precedence — when the caller
            # explicitly provided one, CCD is assumed to be wrong or a
            # generic placeholder (biotite returns one for some codes
            # like "LIG"). Fall through to CCD otherwise.
            if ligand_smiles_dict and resname in ligand_smiles_dict:
                component_smiles = ligand_smiles_dict[resname]
            else:
                reference_code = (ligand_ccd_code_dict or {}).get(
                    resname,
                    resname,
                )
                component_smiles = _get_ccd_smiles(reference_code)
            if component_smiles is None:
                reference_fragments = []
                break
            reference_fragments.append(component_smiles)
        if reference_fragments:
            smiles = ".".join(reference_fragments)
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
        #   - fallback SMILES when the CCD/user-SMILES lookup failed
        resolved_smiles: str | None = None
        stereo_matches: bool | None = None
        try:
            from plinder.data.annotations.cif_utils import atoms_to_rdkit_mol

            # biotite has no chiral tags → stereo assigned from 3D inside helper
            resolved_mol = atoms_to_rdkit_mol(lig_heavy)
            resolved_smiles = str(Chem.MolToSmiles(resolved_mol))
            # Compare resolved 3D stereo with CCD template stereo
            # (works for both single- and multi-residue ligands)
            stereo_matches = _check_stereo_vs_template(
                resolved_mol,
                custom_templates=custom_templates,
                ccd_code_dict=ligand_ccd_code_dict,
            )
        except Exception as e:
            LOG.warning(f"Failed to compute resolved SMILES for {ccd_code}: {e}")
        # Prefer the candidate representing more atoms. Resolved coordinates
        # can join a composite ligand correctly, but may omit unobserved atoms;
        # the complete set of per-component references can be fuller despite
        # not encoding the observed inter-component links. Equal-sized resolved
        # molecules win because they do preserve those links.
        if len(residue_component_ids) > 1:
            smiles = _choose_ligand_smiles_by_heavy_atom_count(
                smiles,
                resolved_smiles,
            )
        elif smiles is None:
            smiles = resolved_smiles
        # Heavy-atom counts are structural (not rdkit descriptors): resolved is
        # just the length of the heavy-atom biotite array (robust even if the
        # rdkit mol build above failed), total comes from the SMILES-keyed
        # descriptor cache, and unresolved is their difference.
        num_resolved_heavy_atoms = lig_heavy.array_length()
        reference_descriptors = _smiles_descriptors(smiles)
        num_heavy_atoms = (
            reference_descriptors["num_heavy_atoms"]
            if reference_descriptors is not None
            else None
        )
        num_unresolved_heavy_atoms = (
            num_heavy_atoms - num_resolved_heavy_atoms
            if num_heavy_atoms and num_resolved_heavy_atoms
            else None
        )
        # Centroid
        centroid = list(lig_atoms.coord.mean(axis=0))
        # BIRD/PRD id straight from the enriched CIF: the mapping key is the PRD
        # code (see cif_utils BIRD parse), a single canonical id for the ligand.
        bird_id = next(iter(ligand_chain.mappings.get("BIRD", {})), "")
        ligand = cls(
            pdb_id=pdb_id,
            biounit_id=biounit_id,
            asym_id=ligand_chain.asym_id,
            instance=ligand_instance,
            ccd_code=ccd_code,
            plip_type=get_chain_type(ligand_chain.chain_type_str),
            bird_id=bird_id,
            centroid=centroid,
            smiles=smiles or "",
            neighboring_residue_threshold=neighboring_residue_threshold,
            neighboring_ligand_threshold=neighboring_ligand_threshold,
            resolved_smiles=resolved_smiles or "",
            resolved_stereo_matches_template=stereo_matches,
            num_heavy_atoms=num_heavy_atoms,
            num_resolved_heavy_atoms=num_resolved_heavy_atoms,
            num_unresolved_heavy_atoms=num_unresolved_heavy_atoms,
            residue_numbers=residue_numbers,
            member_residue_numbers=member_residue_numbers,
        )

        # Find neighboring polymer residues (protein + nucleic acid) without
        # rebuilding or scanning a whole-assembly polymer index.
        neighbor_indices = spatial_index.atom_indices_near(
            lig_coords, ligand.neighboring_residue_threshold
        )
        neighboring_atoms = spatial_index.take_atoms(
            biounit, neighbor_indices, include_bonds=False
        )
        polymer_mask = struc.filter_amino_acids(
            neighboring_atoms
        ) | struc.filter_nucleotides(neighboring_atoms)
        near_prot = neighboring_atoms[polymer_mask]

        for chain_id in np.unique(near_prot.chain_id):
            if chain_id in member_instance_chains:
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
            if c not in member_instance_chains
        }

        # A covalently-linked group is one molecule: bonds *between* member
        # chains are internal, so only count covale bonds from any member to a
        # non-member (receptor) chain.
        neighboring_non_member = neighboring_asym_ids - member_asym_ids
        covalent_linkages: set[str] = set()
        for member_asym in member_asym_ids:
            covalent_linkages |= extract_ligand_links_to_neighbouring_chains(
                all_covalent_dict,
                member_asym,
                neighboring_non_member,
                link_type="covale",
            )
        ligand.covalent_linkages = covalent_linkages
        ligand.is_covalent = len(ligand.covalent_linkages) > 0

        # Find neighboring ligand chains
        near_lig_indices = spatial_index.atom_indices_near(
            lig_coords, ligand.neighboring_ligand_threshold
        )
        near_all = spatial_index.take_atoms(
            biounit, near_lig_indices, include_bonds=False
        )

        ligand.neighboring_ligands = sorted(
            {
                c
                for c in np.unique(near_all.chain_id)
                if c not in member_instance_chains
                and "." in c
                and c.split(".")[1] in ligand_like_chains
            }
        )
        if water_chains is None:
            water_chains = get_water_chain_ids(biounit)
        # Populate interactions and waters from peppr results
        ligand.interactions = peppr_interactions
        ligand.waters = defaultdict(list)
        for w_chain, w_resnum in peppr_waters:
            ligand.waters[w_chain].append(w_resnum)

        # Derive interacting residues from peppr interaction hashes
        for instance_chain, residues in peppr_interactions.items():
            if instance_chain in member_instance_chains:
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
            # obsolete codes are never ingested (PDB remediates), so the old
            # code->canonical synonym map was an identity no-op: use the code as-is.
            ligand.unique_ccd_code = ligand.ccd_code

        return ligand

    @property
    def _members(self) -> dict[str, list[int]]:
        """Residue numbers per member instance-chain (self, incl. covalent).

        Falls back to the single primary chain for ligands built before
        ``member_residue_numbers`` was populated.
        """
        return self.member_residue_numbers or {
            self.instance_chain: self.residue_numbers
        }

    @property
    def _is_multi_residue(self) -> bool:
        """Whether deposited residue membership proves a composite ligand."""
        residue_count = sum(
            len(set(residue_numbers)) for residue_numbers in self._members.values()
        )
        return residue_count > 1 or len(self.ccd_code.split("-")) > 1

    @property
    def member_asym_ids(self) -> list[str]:
        """Asym IDs of every chain this (possibly merged) ligand spans."""
        return sorted({ic.split(".")[-1] for ic in self._members})

    @cached_property
    def selection(self) -> str:
        """[EXCLUDE] Selection string for ligand

        Spans every member instance-chain so covalently-linked ligand chains
        (a macrocycle deposited as several chains) select all of their atoms.
        """

        def _chain_selection(instance_chain: str, resnums: list[int]) -> str:
            selection = f"cname='{instance_chain}'"
            if len(resnums):
                residue_selection = " or ".join(f"rnum={rnum}" for rnum in resnums)
                selection += f"and ({residue_selection})"
            return selection

        members = self._members
        if len(members) == 1:
            ((instance_chain, resnums),) = members.items()
            return _chain_selection(instance_chain, resnums)
        return " or ".join(
            f"({_chain_selection(instance_chain, resnums)})"
            for instance_chain, resnums in members.items()
        )

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
        """[EXCLUDE] Residues in the ligand's binding pocket which includes neighboring and interacting residues."""
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
        if not self.num_heavy_atoms:
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
        """[EXCLUDE] Counter of interactions for a given ligand."""
        interactions_counter: dict[str, dict[int, ty.Counter[str]]] = {}
        for chain in self.interactions:
            interactions_counter[chain] = {}
            for residue in self.interactions[chain]:
                interactions_counter[chain][residue] = Counter(
                    self.interactions[chain][residue]
                )
        return interactions_counter

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
        # Match by CCD code OR by structure (RDKit SMILES), so a ligand that is
        # the same molecule under a different current code is still classified.
        # COFACTOR_SMILES / ARTIFACT_SMILES are None if only the code sets were
        # populated (e.g. a test that patches COFACTORS/ARTIFACTS directly) —
        # fall back to code-only matching then.
        if self.ccd_code in COFACTORS or (
            COFACTOR_SMILES is not None and self.smiles in COFACTOR_SMILES
        ):
            self.is_cofactor = True
        if self.ccd_code in ARTIFACTS or (
            ARTIFACT_SMILES is not None and self.smiles in ARTIFACT_SMILES
        ):
            self.in_artifact_list = True

        if self.is_ion:
            self.is_artifact = False
        elif self.in_artifact_list:
            self.is_artifact = True
        elif lig_has_dummies(self.ccd_code):
            # check for dummy list including composites, too!
            self.is_artifact = True
        elif self._is_multi_residue and any(
            (
                self.is_oligosaccharide,
                self.is_oligonucleotide,
                self.is_oligopeptide,
            )
        ):
            # Small-molecule charge and linker cutoffs do not describe
            # recognized oligomeric ligands.  For example, a short peptide can
            # legitimately exceed the formal-charge cutoff through Lys/Arg.
            self.is_artifact = False
        elif is_excluded_mol(self.smiles):
            self.is_artifact = True
        else:
            self.is_artifact = False

        # reset: artifacts should not count to other ligand class definitions
        if self.is_artifact:
            self.is_ion = False
            self.is_monosaccharide = False
            self.is_oligosaccharide = False
            self.is_mononucleotide = False
            self.is_oligonucleotide = False
            self.is_monopeptide = False
            self.is_oligopeptide = False
            self.is_cofactor = False
            self.is_lipinski = False
            self.is_fragment = False
            self.is_covalent = False

        # Indicator of whether a ligand type is not any recognized class.
        self.is_other = not any(
            [
                self.is_invalid,
                self.is_ion,
                self.is_monosaccharide,
                self.is_oligosaccharide,
                self.is_mononucleotide,
                self.is_oligonucleotide,
                self.is_monopeptide,
                self.is_oligopeptide,
                self.is_artifact,
                self.is_cofactor,
                self.is_lipinski,
                self.is_fragment,
                self.is_covalent,
            ]
        )

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
        List of residues in the format
        ``<instance>.<label_asym_id>_<label_seq_id>_<residue_index>_``
        ``<auth_seq_id>_<insertion_code>``.
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
                    f"{instance_chain}_{residue_number}_{chains[chain].residues[residue_number].index}_{chains[chain].residues[residue_number].auth_number}_{chains[chain].residues[residue_number].insertion_code}"
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
        for field, (description, _) in self.get_descriptions_and_types().items():
            if description_excluded_from_flat_export(description):
                continue
            name = f"ligand_{field}"
            data[name] = getattr(self, field, None)

        # These internal selections are required to reconstruct a system
        # deterministically from the source mmCIF and an annotation row.
        data["ligand_residue_numbers"] = sorted(set(self.residue_numbers))
        # Every instance-chain this ligand spans — a multi-chain covalent ligand
        # spans several — so the source-reconstruction SDF covers the whole
        # molecule, not just the primary chain.
        data["ligand_instance_chains"] = sorted(self._members)
        data["ligand_water_residues"] = sorted(
            f"{chain_id}_{residue_number}"
            for chain_id, residue_numbers in self.waters.items()
            for residue_number in set(residue_numbers)
        )

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
