# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
from __future__ import annotations

import functools
from functools import cached_property
from typing import Any

import biotite.structure as struc
import biotite.structure.io.pdbx as pdbx
from PDBValidation.Validation import PDBValidation
from pydantic import ConfigDict, Field

from plinder.data.utils.annotations.get_ligand_validation import (
    ResidueListValidation,
    ResidueValidation,
    ResidueValidationThresholds,
)
from plinder.data.utils.annotations.utils import DocBaseModel


@functools.cache
def _standard_aa_names() -> set[str]:
    """Standard amino acid 3-letter codes."""
    import biotite.structure.info as info

    return set(info.amino_acid_names())


@functools.cache
def _standard_na_names() -> set[str]:
    """Standard nucleotide 3-letter codes (RNA + DNA)."""
    return {"A", "C", "G", "U", "DA", "DC", "DG", "DT", "DU"}


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


class Residue(DocBaseModel):
    chain: str
    index: int
    number: int
    auth_number: str
    one_letter_code: str
    name: str
    chem_type: str
    validation: ResidueValidation | None = None
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

    Attributes
    ----------
    is_ptm: bool
        Does residue have post translational modification

    Examples
    --------
    TODO: Add example
    """

    @cached_property
    def is_modified(self) -> bool:
        """Is this a modified residue (PTM for protein, modified base for NA)."""
        ct = self.chem_type.lower()
        if "peptide" in ct:
            return self.name not in _standard_aa_names()
        if "rna" in ct or "dna" in ct or "nucleotide" in ct:
            return self.name not in _standard_na_names()
        return False

    @cached_property
    def is_ptm(self) -> bool:
        """Does the residue have a post-translational modification (protein only)."""
        ct = self.chem_type.lower()
        return "peptide" in ct and self.name not in _standard_aa_names()


class Chain(DocBaseModel):
    asym_id: str = Field(description="Chain asymmetric id")
    auth_id: str = Field(description="Chain author id")
    entity_id: str = Field(description="Chain entity id")
    chain_type_str: str = Field(
        description="__Chain type string from CIF entity_poly.type"
    )
    residues: dict[int, Residue] = Field(
        description="__Dictionary of residues in chain with keys as residue number"
    )
    length: int = Field(description="SEQRES length")
    num_unresolved_residues: int = Field(
        description="Number of unresolved residues (SEQRES length - len(residues))"
    )
    mappings: dict[str, dict[str, list[tuple[str, str] | None]]] = Field(
        default_factory=dict,
        description="__Mapping of metadata associated with chain with keys as chain asym id",
    )
    holo: bool = Field(
        default=True, description="__Is the chain part of a holo system or not"
    )
    validation: ResidueListValidation | None = Field(
        default=None,
        description="__Crystal validation information for the residues in the chain",
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
        """
        import biotite.structure as struc
        import biotite.structure.info as info

        # Build residue dict
        residues = {}
        res_starts = struc.get_residue_starts(atoms)
        for idx, start in enumerate(res_starts):
            resnum = int(atoms.res_id[start])
            resname = atoms.res_name[start]
            auth_resnum = str(atoms.res_id[start])
            # One-letter code: try amino acid first, then nucleotide
            olc = "X"
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
            residues[resnum] = Residue(
                chain=asym_id,
                index=idx,
                number=resnum,
                auth_number=auth_resnum,
                one_letter_code=olc,
                name=resname,
                chem_type=chem_type,
            )

        # Get entity_id from _struct_asym
        entity_id = ""
        if "struct_asym" in block:
            sa = block["struct_asym"]
            sa_ids = sa["id"].as_array()
            sa_entities = sa["entity_id"].as_array()
            for i, sa_id in enumerate(sa_ids):
                if sa_id == asym_id:
                    entity_id = sa_entities[i]
                    break

        # Get auth chain ID
        auth_id = ""
        if hasattr(atoms, "auth_asym_id"):
            auth_id = atoms.auth_asym_id[0]
        elif "atom_site" in block:
            atom_site = block["atom_site"]
            if "auth_asym_id" in atom_site:
                asym_arr = atom_site["label_asym_id"].as_array()
                auth_arr = atom_site["auth_asym_id"].as_array()
                for i, a in enumerate(asym_arr):
                    if a == asym_id:
                        auth_id = auth_arr[i]
                        break

        # Get chain type from _entity_poly.type or _entity.type
        chain_type_str = _get_chain_type_from_cif(block, entity_id)

        return cls(
            asym_id=asym_id,
            auth_id=auth_id,
            entity_id=entity_id,
            chain_type_str=chain_type_str,
            residues=residues,
            length=seqres_length,
            num_unresolved_residues=seqres_length - len(residues),
        )

    @cached_property
    def residue_index_to_number(self) -> dict[int, int]:
        """
        __Dictionary of residue index to residue number
        """
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
                self.asym_id, residue, self.entity_id, doc
            )
        validations = []
        for r in self.residues.values():
            if r.validation is not None:
                validations.append(r.validation)
        self.validation = ResidueListValidation.from_residues(
            validations, residue_thresholds
        )
