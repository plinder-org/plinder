"""Prepare coordinate inputs without contact detection or artifact classification."""

from __future__ import annotations

import shutil
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING

import biotite.structure as struc
import biotite.structure.io.pdbx as pdbx
import numpy as np
from biotite.interface.rdkit import to_mol
from rdkit.Chem import Mol

from plinder.core.structure.inputs import (
    StructureInput,
    read_input_structure,
    read_ligand_sdf,
)
from plinder.data.annotations.cif_utils import (
    check_cif_bond_orders,
    check_custom_mmcif_fields,
    enrich_cif_with_ccd_bonds,
    enrich_cif_with_smiles_bonds,
    get_structure_with_altloc,
    read_mmcif_file,
)
from plinder.data.annotations.save_utils import save_cif_file, save_ligands

if TYPE_CHECKING:
    from plinder.core.index.system import PlinderSystem


def reference_ligands(
    reference: PlinderSystem, *, include_all_ligands: bool = False
) -> dict[str, Path]:
    """Select reference poses by release ligand ID, proper ligands by default.

    Uses the reference annotation, not a fresh classification of the prediction.
    A covalent multi-chain ligand remains one complete SDF molecule.
    """
    rows = reference.system
    if not include_all_ligands:
        rows = rows.loc[rows["ligand_is_proper"].eq(True)]
    if rows.empty:
        raise ValueError(f"No selected reference ligands in {reference.system_id}")
    sdfs = reference.ligand_sdfs
    return {
        row.ligand_id: Path(sdfs[row.ligand_instance_chain])
        for row in rows.itertuples()
    }


@dataclass
class PreparedPrediction:
    """Persistent receptor and ligand files, plus a receptor for PoseBusters.

    Ligand keys are supplied SDF filenames or original label asym IDs. Covalently
    connected mmCIF ligand chains form one molecule, keyed by their first asym ID.
    """

    receptor: Path
    ligands: dict[str, Path]
    receptor_molecule: Mol


def _prepare_receptor_and_sdfs(
    model: Path, output_dir: Path, ligand_sdfs: Sequence[str | Path]
) -> PreparedPrediction:
    """Keep supplied poses and load the predicted receptor for PoseBusters."""
    paths = [Path(path).resolve() for path in ligand_sdfs]
    if len({path.name for path in paths}) != len(paths):
        raise ValueError("Supplied ligand SDF filenames must be unique per prediction")
    for path in paths:
        read_ligand_sdf(path)
    atoms, _ = read_input_structure(model)
    receptor = to_mol(atoms[struc.filter_heavy(atoms)])
    output_dir.mkdir(parents=True, exist_ok=True)
    ligands = {}
    for index, path in enumerate(paths):
        destination = output_dir / f"ligand_{index:04d}.sdf"
        if path != destination.resolve():
            shutil.copyfile(path, destination)
        ligands[path.name] = destination.resolve()
    return PreparedPrediction(model.resolve(), ligands, receptor)


def prepare_prediction(
    model: str | Path | StructureInput,
    output_dir: str | Path,
    *,
    ligand_sdfs: Sequence[str | Path] | None = None,
    ligand_smiles: Mapping[str, str] | None = None,
    ligand_ccd_codes: Mapping[str, str] | None = None,
    ligand_chains: Sequence[str] = (),
) -> PreparedPrediction:
    """Prepare a receptor with supplied SDFs, or extract poses from a full mmCIF.

    With ``ligand_sdfs``, ``model`` is the predicted receptor (PDB or mmCIF,
    optionally gzip-compressed). Each SDF contains one ligand in that receptor's
    coordinate frame and supplies its own chemistry. Result IDs use the SDF
    filenames. An empty list represents a prediction with zero ligand poses.
    Chemistry mappings and ``ligand_chains`` apply to complete-mmCIF extraction.

    Complete-mmCIF extraction selects non-water, non-polymer chains as candidate
    ligands, including distant poses, ions and artifacts. Protein and DNA/RNA
    polymers remain receptors. ``ligand_chains`` additionally selects polymer
    chains used as ligands, such as peptides, by their label asym IDs. Without
    entity metadata, chains containing amino-acid or nucleotide residues are
    treated as polymers unless explicitly selected as ligands.
    DNA/RNA inputs need their mmCIF polymer sequence metadata retained.

    CCD/SMILES mappings use component names, e.g. ``{"LIG": "ATP"}``. For
    SMILES, the heavy-atom order must match the CIF atom order; CCD templates
    use atom names or graph matching. Unknown chemistry raises an error with
    the components requiring bond information. This function does not decide
    which ligands to evaluate: select proper ligands on the reference instead.

    Prepared SDFs and extracted receptors remain in ``output_dir``; supplied
    receptors use their original paths. The RDKit receptor retains residue
    metadata for PoseBusters.
    """
    if isinstance(model, StructureInput):
        if ligand_sdfs is not None:
            raise ValueError("Supply ligand_sdfs in the StructureInput only")
        ligand_sdfs = model.ligand_sdfs
        model = model.coordinates
    model, output_dir = Path(model), Path(output_dir)
    if ligand_sdfs is not None:
        return _prepare_receptor_and_sdfs(model, output_dir, ligand_sdfs)
    if model.name.lower().endswith((".pdb", ".pdb.gz")):
        raise ValueError(
            "Ligand preparation requires a complete mmCIF or a receptor with "
            "ligand SDFs. Supply model.sdf or model.ligands/ beside model.pdb."
        )
    cif = read_mmcif_file(model)
    block: pdbx.CIFBlock = next(iter(cif.values()))
    check_custom_mmcif_fields(
        block, source=model, structure_mode="as_is", require_label_ids=True
    )
    overlap = set(ligand_smiles or {}) & set(ligand_ccd_codes or {})
    if overlap:
        raise ValueError(
            f"Provide either SMILES or CCD for each component: {sorted(overlap)}"
        )
    atoms = get_structure_with_altloc(cif)
    atoms = atoms[struc.filter_heavy(atoms) & ~struc.filter_solvent(atoms)]
    chains = set(atoms.chain_id)
    extra_ligands = set(ligand_chains)
    if missing := extra_ligands - chains:
        raise ValueError(f"Ligand chains absent from model: {sorted(missing)}")

    entity_types = {}
    if "entity" in block and "struct_asym" in block:
        entity = block["entity"]
        asym = block["struct_asym"]
        by_entity = dict(zip(entity["id"].as_array(), entity["type"].as_array()))
        entity_types = {
            chain: by_entity.get(entity_id)
            for chain, entity_id in zip(
                asym["id"].as_array(), asym["entity_id"].as_array()
            )
        }
    polymer_atoms = struc.filter_amino_acids(atoms) | struc.filter_nucleotides(atoms)
    receptor_chains = {
        chain
        for chain in chains - extra_ligands
        if entity_types.get(chain) == "polymer"
        or (
            entity_types.get(chain) is None
            and polymer_atoms[atoms.chain_id == chain].any()
        )
    }
    if not receptor_chains:
        raise ValueError(f"No protein or nucleic-acid receptor chains in {model}")
    candidates = chains - receptor_chains
    # The chemistry helpers scan HETATM records. Select by chain membership,
    # not the input record spelling, so ATOM-labelled ligands are checked too.
    site = block["atom_site"]
    site["group_PDB"] = np.where(
        np.isin(site["label_asym_id"].as_array(str), list(candidates)),
        "HETATM",
        site["group_PDB"].as_array(str),
    )
    if ligand_ccd_codes:
        enrich_cif_with_ccd_bonds(cif, dict(ligand_ccd_codes))
    if ligand_smiles:
        enrich_cif_with_smiles_bonds(cif, dict(ligand_smiles))
    check_cif_bond_orders(cif)
    atoms = get_structure_with_altloc(cif, include_bonds=True)
    atoms = atoms[struc.filter_heavy(atoms) & ~struc.filter_solvent(atoms)]
    receptor_atoms = atoms[np.isin(atoms.chain_id, list(receptor_chains))]

    # Join covalently linked ligand chains, never through a receptor or metal.
    neighbors: dict[str, set[str]] = {chain: set() for chain in candidates}
    for first, second, bond_type in atoms.bonds.as_array():
        if bond_type == struc.BondType.COORDINATION:
            continue
        a, b = atoms.chain_id[first], atoms.chain_id[second]
        if a != b and a in candidates and b in candidates:
            neighbors[a].add(b)
            neighbors[b].add(a)
    groups = {}
    unseen = set(candidates)
    while unseen:
        first = min(unseen)
        group, pending = set(), [first]
        while pending:
            chain = pending.pop()
            if chain not in group:
                group.add(chain)
                pending.extend(neighbors[chain] - group)
        unseen -= group
        groups[first] = sorted(group)

    output_dir.mkdir(parents=True, exist_ok=True)
    receptor_path = output_dir / "receptor.cif"
    save_cif_file(
        receptor_atoms,
        "receptor",
        receptor_path,
        source_block=block,
        source_asym_ids={chain: chain for chain in receptor_chains},
    )
    # File names do not interpolate user-supplied chain IDs; keys preserve them.
    names = {chain: f"ligand_{index:04d}" for index, chain in enumerate(groups)}
    save_ligands(
        atoms, {names[chain]: group for chain, group in groups.items()}, output_dir
    )
    return PreparedPrediction(
        receptor=receptor_path,
        ligands={chain: output_dir / f"{name}.sdf" for chain, name in names.items()},
        receptor_molecule=to_mol(receptor_atoms),
    )
