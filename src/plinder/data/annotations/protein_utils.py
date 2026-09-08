# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
from __future__ import annotations

import functools
from collections import Counter
from collections.abc import Iterable, Mapping
from functools import cached_property
from typing import Any, NamedTuple

import biotite.structure as struc
import biotite.structure.io.pdbx as pdbx
import numpy as np
from PDBValidation.Validation import PDBValidation
from pydantic import ConfigDict, Field

from plinder.core.utils.log import setup_logger
from plinder.data.annotations.get_ligand_validation import (
    ResidueListValidation,
    ResidueValidation,
    ResidueValidationThresholds,
)
from plinder.data.annotations.utils import DocBaseModel


@functools.cache
def _standard_aa_names() -> set[str]:
    """Standard amino acid 3-letter codes."""
    import biotite.structure.info as info

    return set(info.amino_acid_names())


@functools.cache
def _standard_na_names() -> set[str]:
    """Standard nucleotide 3-letter codes (RNA + DNA)."""
    return {"A", "C", "G", "U", "DA", "DC", "DG", "DT", "DU"}


@functools.cache
def _ccd_parent_components() -> dict[str, str]:
    """``comp_id -> mon_nstd_parent_comp_id`` from biotite's bundled CCD."""
    from biotite.structure.info import ccd as bundled_ccd

    chem_comp = bundled_ccd.get_ccd()["chem_comp"]
    return dict(
        zip(
            chem_comp["id"].as_array(str),
            chem_comp["mon_nstd_parent_comp_id"].as_array(str),
        )
    )


LOG = setup_logger(__name__)

_CANONICAL_MONOMERS = frozenset(
    {
        "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
        "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
        "UNK", "A", "C", "G", "U", "I", "N", "DA", "DC", "DG", "DT", "DU", "DI", "DN",
    }
)  # fmt: skip


def get_modified_residues(
    block: pdbx.CIFBlock,
    residue_author_ids: Mapping[str, Mapping[int, tuple[str, str]]] | None = None,
) -> dict[str, list[str]]:
    """Non-canonical SEQRES monomers per label asym.

    Entries are the :func:`~plinder.data.annotations.cif_utils.residue_address`
    ``{auth_seq}:{comp_id}:{asym}:{label_seq}`` followed by ``>{parent}`` and
    `` (details)`` from ``pdbx_struct_mod_residue`` when present. Built from
    ``entity_poly_seq`` (``comp_id`` is its ``mon_id``), so unresolved positions
    are included with ``?`` as ``auth_seq``. Parent from
    ``pdbx_struct_mod_residue``, else the CCD, else ``?``.
    """
    from plinder.data.annotations.cif_utils import _iter_category_rows, residue_address

    if residue_author_ids is None:
        _, residue_author_ids = get_atom_site_author_ids(block)
    modified_by_entity: dict[str, list[tuple[int, str]]] = {}
    for row in _iter_category_rows(
        block, "entity_poly_seq", ["entity_id", "num", "mon_id"]
    ):
        if row["mon_id"] in _CANONICAL_MONOMERS:
            continue
        try:
            number = int(row["num"])
        except ValueError:
            continue
        modified_by_entity.setdefault(row["entity_id"], []).append(
            (number, row["mon_id"])
        )
    if not modified_by_entity:
        return {}
    missing = {"", ".", "?"}
    deposited: dict[tuple[str, int, str], tuple[str, str]] = {}
    for row in _iter_category_rows(
        block,
        "pdbx_struct_mod_residue",
        ["label_asym_id", "label_seq_id", "label_comp_id", "parent_comp_id", "details"],
    ):
        try:
            key = (row["label_asym_id"], int(row["label_seq_id"]), row["label_comp_id"])
        except ValueError:
            continue
        deposited[key] = (row["parent_comp_id"], row["details"])
    ccd_parents = _ccd_parent_components()
    modified: dict[str, list[str]] = {}
    for row in _iter_category_rows(block, "struct_asym", ["id", "entity_id"]):
        positions = modified_by_entity.get(row["entity_id"])
        if not positions:
            continue
        asym_id = row["id"]
        author_ids = residue_author_ids.get(asym_id, {})
        entries = []
        for number, mon_id in sorted(positions):
            parent, details = deposited.get((asym_id, number, mon_id), ("", ""))
            if parent in missing:
                parent = ccd_parents.get(mon_id, "?")
            if parent in missing:
                parent = "?"
            auth_seq, insertion = author_ids.get(number, ("?", "."))
            if insertion not in missing:
                auth_seq = f"{auth_seq}{insertion}"
            entry = f"{residue_address(auth_seq, mon_id, asym_id, number)}>{parent}"
            if details not in missing:
                entry = f"{entry} ({details})"
            entries.append(entry)
        modified[asym_id] = entries
    return modified


class UnobservedAtom(NamedTuple):
    """A ``pdbx_unobs_or_zero_occ_atoms`` row in ``residue_address`` fields."""

    comp_id: str
    auth_seq: str  # auth_seq_id with insertion code
    label_seq: str  # label_seq_id; "." for non-polymer and branched residues
    atom_name: str
    zero_occupancy: bool  # occupancy_flag 0: modelled at zero occupancy, not absent


def get_unobserved_atoms(
    block: pdbx.CIFBlock,
) -> tuple[dict[str, dict[int, list[str]]], dict[str, list[UnobservedAtom]]] | None:
    """Model-1 rows of ``pdbx_unobs_or_zero_occ_atoms``, or ``None`` when absent.

    Returns ``(asym -> label_seq_id -> atom names)`` for polymer residues and
    ``(asym -> [UnobservedAtom])`` for every row.
    """
    from plinder.data.annotations.cif_utils import _iter_category_rows

    category = "pdbx_unobs_or_zero_occ_atoms"
    if category not in block:
        return None
    columns = [
        "PDB_model_num",
        "polymer_flag",
        "label_asym_id",
        "label_comp_id",
        "label_seq_id",
        "auth_seq_id",
        "PDB_ins_code",
        "label_atom_id",
        "occupancy_flag",
    ]
    absent = [column for column in columns if column not in block[category]]
    if absent:
        LOG.warning(f"{category} lacks {absent}; ignoring the category")
        return None
    missing = {"", ".", "?"}
    by_residue: dict[str, dict[int, list[str]]] = {}
    by_chain: dict[str, list[UnobservedAtom]] = {}
    for row in _iter_category_rows(block, category, columns):
        if row["PDB_model_num"] not in {"1", "?", "."}:
            continue
        asym_id = row["label_asym_id"]
        auth_seq = row["auth_seq_id"]
        if row["PDB_ins_code"] not in missing:
            auth_seq = f"{auth_seq}{row['PDB_ins_code']}"
        is_polymer = row["polymer_flag"] == "Y"
        by_chain.setdefault(asym_id, []).append(
            UnobservedAtom(
                row["label_comp_id"],
                auth_seq,
                row["label_seq_id"] if is_polymer else ".",
                row["label_atom_id"],
                row["occupancy_flag"] == "0",
            )
        )
        try:
            number = int(row["label_seq_id"])
        except ValueError:
            continue
        by_residue.setdefault(asym_id, {}).setdefault(number, []).append(
            row["label_atom_id"]
        )
    return by_residue, by_chain


def _get_chain_type_from_cif(block: pdbx.CIFBlock, entity_id: str) -> str:
    """Get chain type string from CIF entity/entity_poly categories."""
    # Try _entity_poly.type first
    if "entity_poly" in block:
        ep = block["entity_poly"]
        ep_ids = ep["entity_id"].as_array()
        ep_types = ep["type"].as_array()
        for i, eid in enumerate(ep_ids):
            if eid == entity_id:
                return str(ep_types[i])
    # Fall back to _entity.type
    if "entity" in block:
        ent = block["entity"]
        ent_ids = ent["id"].as_array()
        ent_types = ent["type"].as_array()
        for i, eid in enumerate(ent_ids):
            if eid == entity_id:
                return str(ent_types[i])
    return "unknown"


def get_atom_site_author_ids(
    block: pdbx.CIFBlock,
) -> tuple[dict[str, str], dict[str, dict[int, tuple[str, str]]]]:
    """Return author chain and residue IDs keyed by label-asym coordinates.

    The residue mapping contains ``auth_seq_id`` and ``pdbx_PDB_ins_code`` for
    each resolved ``(label_asym_id, label_seq_id)`` pair.  Custom mmCIFs may
    omit author fields, in which case the corresponding label IDs are used.
    """
    if "atom_site" not in block:
        return {}, {}
    atom_site = block["atom_site"]
    if not {"label_asym_id", "label_seq_id"}.issubset(atom_site):
        return {}, {}

    label_asym_ids = atom_site["label_asym_id"].as_array()
    label_seq_ids = atom_site["label_seq_id"].as_array()
    auth_asym_ids = (
        atom_site["auth_asym_id"].as_array()
        if "auth_asym_id" in atom_site
        else label_asym_ids
    )
    auth_seq_ids = (
        atom_site["auth_seq_id"].as_array()
        if "auth_seq_id" in atom_site
        else label_seq_ids
    )
    insertion_codes = (
        atom_site["pdbx_PDB_ins_code"].as_array()
        if "pdbx_PDB_ins_code" in atom_site
        else np.full(len(label_asym_ids), ".", dtype=object)
    )

    chain_auth_ids: dict[str, str] = {}
    residue_author_ids: dict[str, dict[int, tuple[str, str]]] = {}
    missing = {"", ".", "?"}
    for label_asym, label_seq, auth_asym, auth_seq, insertion_code in zip(
        label_asym_ids,
        label_seq_ids,
        auth_asym_ids,
        auth_seq_ids,
        insertion_codes,
    ):
        asym_id = str(label_asym)
        author_chain = str(auth_asym)
        if author_chain in missing:
            author_chain = asym_id
        chain_auth_ids.setdefault(asym_id, author_chain)
        try:
            residue_number = int(str(label_seq))
        except ValueError:
            # Non-polymer atoms use missing label_seq_id values and are not
            # part of a receptor residue mapping.
            continue
        author_number = str(auth_seq)
        if author_number in missing:
            author_number = str(residue_number)
        insertion = str(insertion_code)
        if insertion in missing:
            insertion = "."
        residue_author_ids.setdefault(asym_id, {}).setdefault(
            residue_number, (author_number, insertion)
        )
    return chain_auth_ids, residue_author_ids


def get_seqres_from_cif(block: pdbx.CIFBlock) -> dict[str, str]:
    """Extract SEQRES (one-letter sequences) per chain from CIF."""
    seqres: dict[str, str] = {}
    if "entity_poly" not in block:
        return seqres
    ep = block["entity_poly"]
    if "pdbx_strand_id" not in ep or "pdbx_seq_one_letter_code_can" not in ep:
        return seqres
    strand_ids = ep["pdbx_strand_id"].as_array()
    sequences = ep["pdbx_seq_one_letter_code_can"].as_array()
    for strands, seq in zip(strand_ids, sequences):
        # Clean up sequence (remove newlines, semicolons)
        clean_seq = seq.replace("\n", "").replace(";", "").strip()
        for chain_id in strands.split(","):
            seqres[chain_id.strip()] = clean_seq
    return seqres


def _is_polypeptide(chain_type_str: str) -> bool:
    return "polypeptide" in chain_type_str.lower()


def _is_polynucleotide(chain_type_str: str) -> bool:
    return (
        "polyribonucleotide" in chain_type_str.lower()
        or "polydeoxyribonucleotide" in chain_type_str.lower()
    )


def get_receptor_type(chain_types: Iterable[str]) -> str:
    """Return a deterministic receptor-composition label for CIF chain types."""
    components: set[str] = set()
    for chain_type in chain_types:
        normalized = chain_type.lower()
        if "polypeptide" in normalized:
            components.add("protein")
        if "polydeoxyribonucleotide" in normalized:
            components.add("dna")
        if "polyribonucleotide" in normalized:
            components.add("rna")
        if not any(
            marker in normalized
            for marker in (
                "polypeptide",
                "polydeoxyribonucleotide",
                "polyribonucleotide",
            )
        ):
            components.add("other")
    order = ("protein", "dna", "rna", "other")
    return "+".join(component for component in order if component in components)


def _is_polysaccharide(chain_type_str: str) -> bool:
    return (
        "polysaccharide" in chain_type_str.lower()
        or "oligosaccharide" in chain_type_str.lower()
        or "branched" in chain_type_str.lower()
    )


def _is_water(chain_type_str: str) -> bool:
    return "water" in chain_type_str.lower()


def _is_polymer(chain_type_str: str) -> bool:
    return "poly" in chain_type_str.lower()


def sequences_match_core(seq_a: str, seq_b: str, min_coverage: float = 0.9) -> bool:
    """Check that two sequences share an identical core (no internal mutations).

    Allows terminal overhangs (N/C-term tags, signal peptides, construct
    boundaries) but rejects any substitution in the aligned region.

    Uses local alignment to find the best-scoring overlap, then verifies
    that every aligned position is identical and the alignment covers at
    least *min_coverage* of the shorter sequence.

    Parameters
    ----------
    seq_a, seq_b : str
        Protein sequences to compare.
    min_coverage : float
        Minimum fraction of the shorter sequence that must be aligned.

    Returns
    -------
    bool
        True if the core overlap is 100% identical and coverage is sufficient.
    """
    from biotite.sequence import ProteinSequence
    from biotite.sequence.align import SubstitutionMatrix, align_optimal

    if not seq_a or not seq_b:
        return False
    try:
        s1 = ProteinSequence(seq_a)
        s2 = ProteinSequence(seq_b)
    except Exception:
        return False
    matrix = SubstitutionMatrix.std_protein_matrix()
    alignments = align_optimal(s1, s2, matrix, local=True)
    if not alignments:
        return False
    trace = alignments[0].trace
    n_aligned = 0
    n_identical = 0
    for i, j in trace:
        if i != -1 and j != -1:
            n_aligned += 1
            if s1[i] == s2[j]:
                n_identical += 1
    min_len = min(len(s1), len(s2))
    return n_aligned == n_identical and n_aligned >= min_coverage * min_len


def detect_ligand_chains(
    entry: Any,
    min_polymer_size: int = 12,
) -> dict[str, str]:
    """Detect which chains are ligands vs receptor polymers.

    A polymer chain (protein, NA, saccharide) with >= min_polymer_size
    residues is receptor.  Everything else — non-polymers, short
    polymers, and BIRD-annotated chains — is a ligand.

    Default threshold of 12 is the minimum length for meaningful
    sequence searches (MMseqs2/Foldseek).
    """
    ligand_chains = dict()
    for chain_name, chain in entry.chains.items():
        ct = chain.chain_type_str
        if _is_water(ct):
            continue

        chain_length = len(chain.residues)
        bird_id = list(chain.mappings.get("BIRD", {"": None}))[0]

        # BIRD-annotated short chains are ligands irrespective of polymer type or length
        if bird_id:
            ligand_chains[chain_name] = ct
        # Polymers >= threshold are receptor
        elif _is_polymer(ct) and chain_length >= min_polymer_size:
            continue
        # Everything else is ligand
        else:
            ligand_chains[chain_name] = ct
    return ligand_chains


def detect_ligand_chains_from_cif(
    block: pdbx.CIFBlock,
    min_polymer_size: int = 12,
) -> dict[str, str] | None:
    """Preflight ligand-chain detection without building a bonded AtomArray.

    The residue count mirrors :func:`detect_ligand_chains`: consecutive
    ``label_seq_id`` values in model 1 correspond to the residue starts used
    by :class:`Chain`. ``None`` means the CIF lacks the columns needed for a
    reliable decision and the caller must fall back to full structure loading.
    """
    if "struct_asym" not in block or "atom_site" not in block:
        return None
    struct_asym = block["struct_asym"]
    atom_site = block["atom_site"]
    if not {"id", "entity_id"}.issubset(struct_asym) or not {
        "label_asym_id",
        "label_seq_id",
    }.issubset(atom_site):
        return None

    atom_asym_ids = np.asarray(atom_site["label_asym_id"].as_array(), dtype=str)
    atom_seq_ids = np.asarray(atom_site["label_seq_id"].as_array(), dtype=str)
    if "pdbx_PDB_model_num" in atom_site:
        model_numbers = np.asarray(
            atom_site["pdbx_PDB_model_num"].as_array(), dtype=str
        )
        model_mask = model_numbers == "1"
        atom_asym_ids = atom_asym_ids[model_mask]
        atom_seq_ids = atom_seq_ids[model_mask]

    valid_residue = ~np.isin(atom_seq_ids, (".", "?", ""))
    residue_asym_ids = atom_asym_ids[valid_residue]
    residue_seq_ids = atom_seq_ids[valid_residue]
    residue_counts: Counter[str] = Counter()
    if len(residue_asym_ids):
        residue_starts = np.ones(len(residue_asym_ids), dtype=bool)
        residue_starts[1:] = (residue_asym_ids[1:] != residue_asym_ids[:-1]) | (
            residue_seq_ids[1:] != residue_seq_ids[:-1]
        )
        residue_counts.update(residue_asym_ids[residue_starts])

    bird_asym_ids: set[str] = set()
    if "pdbx_molecule" in block and "asym_id" in block["pdbx_molecule"]:
        bird_asym_ids.update(
            str(value) for value in block["pdbx_molecule"]["asym_id"].as_array()
        )

    observed_asym_ids = set(atom_asym_ids)
    ligand_chains: dict[str, str] = {}
    for asym_id, entity_id in zip(
        struct_asym["id"].as_array(),
        struct_asym["entity_id"].as_array(),
    ):
        asym_id = str(asym_id)
        if asym_id not in observed_asym_ids:
            continue
        chain_type = _get_chain_type_from_cif(block, str(entity_id))
        if _is_water(chain_type):
            continue
        if asym_id in bird_asym_ids:
            ligand_chains[asym_id] = chain_type
        elif _is_polymer(chain_type) and residue_counts[asym_id] >= min_polymer_size:
            continue
        else:
            ligand_chains[asym_id] = chain_type
    return ligand_chains


class Residue(DocBaseModel):
    chain: str
    index: int
    number: int
    auth_number: str
    insertion_code: str = "."
    one_letter_code: str
    name: str
    chem_type: str
    validation: ResidueValidation | None = None
    selected_altcode: str = Field(
        default=".",
        exclude=True,
        description="[EXCLUDE] Deposited alternate conformer selected for this residue",
    )
    unresolved_atom_names: list[str] = Field(
        default_factory=list,
        description="[EXCLUDE] label_atom_id of heavy atoms missing from the model (wwPDB unobserved-atom records, else CCD template minus resolved atoms)",
    )
    """Single residue in a polymer chain.

    Parameters
    ----------
    chain : str
        chain name
    index : int
        residue index
    number : int
        residue number
    one_letter_code : str
        residue one-letter code
    name : str
        residue name
    chem_type: str
        residue chemical type

    Examples
    --------
    TODO: Add example
    """


class Chain(DocBaseModel):
    asym_id: str = Field(description="Chain asymmetric id")
    auth_id: str = Field(description="Chain author id")
    entity_id: str = Field(description="Chain entity id")
    chain_type_str: str = Field(
        description="[EXCLUDE] Chain type string from CIF entity_poly.type"
    )
    residues: dict[int, Residue] = Field(
        description="[EXCLUDE] Dictionary of residues in chain with keys as residue number"
    )
    length: int = Field(description="SEQRES length")
    num_unresolved_residues: int = Field(
        description="Number of unresolved residues (SEQRES length - len(residues))"
    )
    mappings: dict[str, dict[str, list[tuple[str, str] | None]]] = Field(
        default_factory=dict,
        description="[EXCLUDE] Mapping of metadata associated with chain with keys as chain asym id",
    )
    holo: bool = Field(
        default=True, description="[EXCLUDE] Is the chain part of a holo system or not"
    )
    validation: ResidueListValidation | None = Field(
        default=None,
        description="[EXCLUDE] Crystal validation information for the residues in the chain",
    )
    modified_residues: list[str] = Field(
        default_factory=list,
        description="[EXCLUDE] Non-canonical SEQRES monomers as {auth_seq}:{comp_id}:{asym}:{label_seq}>{parent} (details); residue address as in ligand_covalent_linkages",
    )

    # Allow arbitrary types for cached properties
    model_config = ConfigDict(
        arbitrary_types_allowed=True,
    )

    @classmethod
    def from_cif_data(
        cls,
        asym_id: str,
        block: "pdbx.CIFBlock",
        atoms: "struc.AtomArray",
        seqres_length: int,
        *,
        entity_id: str | None = None,
        auth_id: str | None = None,
        chain_type_str: str | None = None,
        residue_author_ids: Mapping[int, tuple[str, str]] | None = None,
        modified_residues: list[str] | None = None,
        unresolved_atoms: Mapping[int, list[str]] | None = None,
    ) -> "Chain":
        """Create Chain from biotite CIF data.

        Parameters
        ----------
        asym_id : str
            Chain asymmetric ID.
        block : pdbx.CIFBlock
            CIF data block for metadata lookup.
        atoms : AtomArray
            Atoms belonging to this chain.
        seqres_length : int
            SEQRES length.
        entity_id, auth_id, chain_type_str : str, optional
            Pre-indexed CIF metadata.  Supplying these avoids repeatedly
            scanning large ``struct_asym`` and entity categories when an
            entry contains thousands of chains.
        residue_author_ids : mapping, optional
            Pre-indexed ``label_seq_id -> (auth_seq_id, insertion_code)``
            mapping for this label-asym chain.
        modified_residues : list[str], optional
            Pre-indexed :func:`get_modified_residues` entries for this chain.
        unresolved_atoms : mapping, optional
            ``label_seq_id -> atom names`` from the wwPDB unobserved-atom records;
            ``None`` derives them from the CCD template of each residue instead.
        """
        import biotite.structure as struc
        import biotite.structure.info as info

        from plinder.core.structure.ccd_template import unresolved_atoms_from_template

        if auth_id is None or residue_author_ids is None:
            chain_auth_ids, author_residues_by_chain = get_atom_site_author_ids(block)
            if auth_id is None:
                auth_id = chain_auth_ids.get(asym_id)
            if residue_author_ids is None:
                residue_author_ids = author_residues_by_chain.get(asym_id, {})

        # Build residue dict
        residues = {}
        res_starts = struc.get_residue_starts(atoms, add_exclusive_stop=True)
        for idx, (start, stop) in enumerate(zip(res_starts[:-1], res_starts[1:])):
            resnum = int(atoms.res_id[start])
            resname = atoms.res_name[start]
            auth_resnum, insertion_code = residue_author_ids.get(
                resnum, (str(resnum), ".")
            )
            selected_altcode = "."
            if hasattr(atoms, "selected_altloc_id"):
                selected_altcode = next(
                    (
                        str(altcode)
                        for altcode in atoms.selected_altloc_id[start:stop]
                        if str(altcode) not in {".", "?", " ", ""}
                    ),
                    ".",
                )
            # One-letter code: try amino acid first, then nucleotide
            olc = "X"
            if chain_type_str is None or "non-polymer" not in chain_type_str.lower():
                try:
                    olc_aa = info.one_letter_code(resname)
                    if olc_aa is not None:
                        olc = olc_aa
                except Exception:
                    pass
            if olc == "X" and resname in _standard_na_names():
                # Map standard nucleotides to their base letter
                olc = resname[-1] if len(resname) <= 2 else resname[1]
            # Determine chem_type from residue name
            if resname in _standard_aa_names():
                chem_type = "Peptide Linking"
            elif resname in _standard_na_names():
                chem_type = (
                    "RNA Linking" if resname in {"A", "C", "G", "U"} else "DNA Linking"
                )
            else:
                chem_type = "Non-Polymer"
            if unresolved_atoms is not None:
                unresolved = list(unresolved_atoms.get(resnum, []))
            else:
                unresolved = (
                    unresolved_atoms_from_template(
                        resname, atoms.atom_name[start:stop].tolist()
                    )
                    or []
                )
            residues[resnum] = Residue(
                chain=asym_id,
                index=idx,
                number=resnum,
                auth_number=auth_resnum,
                insertion_code=insertion_code,
                one_letter_code=olc,
                name=resname,
                chem_type=chem_type,
                selected_altcode=selected_altcode,
                unresolved_atom_names=unresolved,
            )

        # Get entity_id from _struct_asym
        if entity_id is None and "struct_asym" in block:
            sa = block["struct_asym"]
            sa_ids = sa["id"].as_array()
            sa_entities = sa["entity_id"].as_array()
            for i, sa_id in enumerate(sa_ids):
                if sa_id == asym_id:
                    entity_id = str(sa_entities[i])
                    break
        if entity_id is None:
            entity_id = ""

        # Get auth chain ID
        if auth_id is None and hasattr(atoms, "auth_asym_id"):
            auth_id = str(atoms.auth_asym_id[0])
        elif auth_id is None and "atom_site" in block:
            atom_site = block["atom_site"]
            if "auth_asym_id" in atom_site:
                asym_arr = atom_site["label_asym_id"].as_array()
                auth_arr = atom_site["auth_asym_id"].as_array()
                for i, a in enumerate(asym_arr):
                    if a == asym_id:
                        auth_id = str(auth_arr[i])
                        break
        if auth_id is None:
            auth_id = ""

        # Get chain type from _entity_poly.type or _entity.type
        if chain_type_str is None:
            chain_type_str = _get_chain_type_from_cif(block, entity_id)

        return cls(
            asym_id=asym_id,
            auth_id=auth_id,
            entity_id=entity_id,
            chain_type_str=chain_type_str,
            residues=residues,
            length=seqres_length,
            # -1 sentinel (collation flags <0) when SEQRES is missing or
            # shorter than the resolved residues — surface it
            num_unresolved_residues=max(-1, seqres_length - len(residues)),
            modified_residues=list(modified_residues or []),
        )

    @cached_property
    def residue_index_to_number(self) -> dict[int, int]:
        """[EXCLUDE] Dictionary of residue index to residue number"""
        return {self.residues[r].index: r for r in self.residues}

    def format(self, instance: int) -> dict[str, Any]:
        """
        Format chain as a dictionary
        """
        data: dict[str, Any] = {
            "asym_id": f"{instance}.{self.asym_id}",
            "auth_id": self.auth_id,
            "entity_id": self.entity_id,
            "length": self.length,
            "num_unresolved_residues": self.num_unresolved_residues,
        }
        # data.update(self.mappings)
        if self.validation is not None:
            data.update(self.validation.format())
        return data

    def set_validation(
        self, doc: PDBValidation, residue_thresholds: ResidueValidationThresholds
    ) -> None:
        for residue in self.residues:
            self.residues[residue].validation = ResidueValidation.from_residue(
                self.asym_id,
                residue,
                self.entity_id,
                doc,
                preferred_altcode=self.residues[residue].selected_altcode,
            )
        validations = []
        for r in self.residues.values():
            if r.validation is not None:
                validations.append(r.validation)
        self.validation = ResidueListValidation.from_residues(
            validations, residue_thresholds
        )
