# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
"""Lightweight entry/system/chain views loaded from the index parquet.

Consumed by similarity scoring.
"""
from __future__ import annotations

from collections import Counter, defaultdict
from collections.abc import Iterable as IterableABC
from dataclasses import dataclass, field
from functools import cached_property
from pathlib import Path
from typing import Any, Iterable

import pandas as pd

from plinder.core.scores.index import query_index
from plinder.core.scores.query import FILTER
from plinder.core.utils.log import setup_logger

LOG = setup_logger(__name__)


@dataclass(frozen=True)
class ChainView:
    asym_id: str  # label_asym_id (e.g. "A")
    auth_id: str  # auth_asym_id (e.g. "A")
    length: int  # SEQRES length
    entity_id: str = ""
    chain_type: str = "polypeptide(L)"
    holo: bool = True
    uniprot_ids: tuple[str, ...] = ()

    @property
    def is_polypeptide(self) -> bool:
        return "polypeptide" in self.chain_type.lower()

    @property
    def receptor_type(self) -> str:
        """Receptor polymer category derived from the CIF chain type."""
        normalized = self.chain_type.lower()
        components = []
        if "polypeptide" in normalized:
            components.append("protein")
        if "polydeoxyribonucleotide" in normalized:
            components.append("dna")
        if "polyribonucleotide" in normalized:
            components.append("rna")
        return "+".join(components) or "other"


@dataclass
class LigandView:
    """Ligand-level annotations needed by pairwise similarity scoring."""

    id: str
    pdb_id: str
    system_id: str
    instance_chain: str
    asym_id: str
    is_proper: bool
    protein_chains_asym_id: list[str]  # protein-only receptor instance chains
    num_pocket_residues: int
    num_interactions: int
    num_unique_interactions: int
    is_3d_score_able: bool = True
    # receptor instance_chain -> {residue_number: residue_index}
    pocket_residue_number_to_index: dict[str, dict[int, int]] = field(
        default_factory=dict
    )
    # receptor instance_chain -> {residue_number: Counter[interaction_type]}
    interactions_counter: dict[str, dict[int, Counter[str]]] = field(
        default_factory=dict
    )


@dataclass
class SystemView:
    id: str
    pdb_id: str
    system_type: str
    protein_chains_asym_id: list[str]  # protein-only, e.g. ["1.A", "1.B"]
    proper_num_pocket_residues: int
    proper_num_interactions: int
    proper_num_unique_interactions: int
    receptor_type: str = "protein"
    # instance_chain -> {residue_number: residue_index}
    pocket_residue_number_to_index: dict[str, dict[int, int]] = field(
        default_factory=dict
    )
    # instance_chain -> {residue_number: Counter[interaction_type]}
    interactions_counter: dict[str, dict[int, Counter[str]]] = field(
        default_factory=dict
    )
    ligands: dict[str, LigandView] = field(default_factory=dict)

    @cached_property
    def pocket_residue_index_to_number(self) -> dict[str, dict[int, int]]:
        return {
            chain: {idx: num for num, idx in mapping.items()}
            for chain, mapping in self.pocket_residue_number_to_index.items()
        }


@dataclass(frozen=True)
class InterfaceView:
    """One unordered protein-chain interface in a biological assembly."""

    id: str
    pdb_id: str
    biounit_id: str
    chain_1: str
    chain_2: str
    chain_1_residue_number_to_index: dict[int, int]
    chain_2_residue_number_to_index: dict[int, int]
    num_contact_residue_pairs: int

    @property
    def chains(self) -> tuple[str, str]:
        return self.chain_1, self.chain_2

    @property
    def residue_number_to_index(self) -> dict[str, dict[int, int]]:
        return {
            self.chain_1: self.chain_1_residue_number_to_index,
            self.chain_2: self.chain_2_residue_number_to_index,
        }


@dataclass
class EntryView:
    pdb_id: str
    chains: dict[str, ChainView]  # by asym_id
    systems: dict[str, SystemView]  # by system_id
    author_to_asym: dict[str, str]
    interfaces: dict[str, InterfaceView] = field(default_factory=dict)

    @cached_property
    def pocket_index_to_number_per_chain(self) -> dict[str, dict[int, int]]:
        """Union of pocket (residue_index → residue_number) maps across every
        system / instance that touches each chain. Used by Scorer's map_row
        to translate foldseek's 0-based residue indices into PDB residue
        numbers so downstream scoring is source-agnostic."""
        result: dict[str, dict[int, int]] = {}
        for system in self.systems.values():
            for instance_chain, n2i in system.pocket_residue_number_to_index.items():
                asym = instance_chain.split(".", 1)[1]
                bucket = result.setdefault(asym, {})
                for num, idx in n2i.items():
                    bucket[idx] = num
        return result

    @cached_property
    def selected_index_to_number_per_chain(self) -> dict[str, dict[int, int]]:
        """Union of ligand-pocket and protein-interface residue mappings.

        Compact mapped alignments retain only these selected query positions.
        This union supports exact ligand pocket and interface reconstruction
        from one release artifact without storing full residue alignments.
        """
        result = {
            asym_id: dict(index_to_number)
            for asym_id, index_to_number in self.pocket_index_to_number_per_chain.items()
        }
        for interface in self.interfaces.values():
            for (
                instance_chain,
                number_to_index,
            ) in interface.residue_number_to_index.items():
                asym_id = instance_chain.split(".", maxsplit=1)[-1]
                bucket = result.setdefault(asym_id, {})
                for number, index in number_to_index.items():
                    previous = bucket.setdefault(index, number)
                    if previous != number:
                        raise ValueError(
                            f"conflicting residue number for {self.pdb_id} "
                            f"chain {asym_id} index {index}: {previous} != {number}"
                        )
        return result

    def chains_for_alignment(self, chain_type: str, aln_type: str) -> list[str]:
        if chain_type not in {"apo", "holo", "pred"}:
            raise ValueError(f"unknown chain_type={chain_type!r}")
        if aln_type not in {"foldseek", "mmseqs"}:
            raise ValueError(f"unknown aln_type={aln_type!r}")

        if chain_type == "holo":
            receptor_asym_ids = {
                instance_chain.split(".", 1)[1]
                for system in self.systems.values()
                if system.system_type == "holo"
                for instance_chain in system.protein_chains_asym_id
            }
            receptor_asym_ids.update(
                instance_chain.split(".", maxsplit=1)[-1]
                for interface in self.interfaces.values()
                for instance_chain in interface.chains
            )
            chains = sorted(
                self.chains[asym].auth_id
                for asym in receptor_asym_ids
                if asym in self.chains and self.chains[asym].is_polypeptide
            )
        elif chain_type == "apo":
            holo_entities = {
                chain.entity_id for chain in self.chains.values() if chain.holo
            }
            chains = sorted(
                chain.auth_id
                for chain in self.chains.values()
                if not chain.holo
                and chain.entity_id not in holo_entities
                and chain.is_polypeptide
            )
        else:
            uniprot_ids = sorted(
                {
                    uniprot_id
                    for chain in self.chains.values()
                    if chain.holo and chain.is_polypeptide
                    for uniprot_id in chain.uniprot_ids
                }
            )
            if aln_type == "foldseek":
                return [f"AF-{uniprot_id}-F1-model_v4_A" for uniprot_id in uniprot_ids]
            return uniprot_ids

        if aln_type == "foldseek":
            return [f"pdb_0000{self.pdb_id}_xyz-enrich_{c}" for c in chains]
        return [f"{self.pdb_id}_{c}" for c in chains]


def _parse_neighboring_residue(s: str) -> tuple[str, int, int]:
    """Parse ``{instance_chain}_{res_number}_{res_index}_{auth_number}`` strings."""
    parts = s.split("_")
    # instance_chain has a dot ("1.A") so it doesn't collide with the _ separator
    return parts[0], int(parts[1]), int(parts[2])


def _parse_interaction(s: str) -> tuple[str, int, str]:
    """Parse ``{instance_chain}_{res_number}_{interaction_type}`` strings.

    interaction_type contains underscores (e.g. ``type:hydrogen_bonds__...``),
    so use maxsplit=2.
    """
    inst, rnum, itype = s.split("_", 2)
    return inst, int(rnum), itype


def _as_list(value: Any) -> list[Any]:
    """Convert an array-valued cell (numpy array / None / NaN) to a list."""
    if value is None:
        return []
    if isinstance(value, str):
        return [value]
    if isinstance(value, IterableABC):
        return list(value)
    return []


def _entry_chains_from_rows(
    entry_rows: pd.DataFrame,
    chain_rows: pd.DataFrame | None = None,
) -> tuple[dict[str, ChainView], dict[str, str]]:
    """Build chain views from the chain table or holo fallback."""
    chains: dict[str, ChainView] = {}
    author_to_asym: dict[str, str] = {}
    if chain_rows is not None:
        for row in chain_rows.itertuples(index=False):
            chain = ChainView(
                asym_id=str(row.chain_asym_id),
                auth_id=str(row.chain_auth_id),
                entity_id=str(row.chain_entity_id),
                chain_type=str(row.chain_type),
                length=int(row.chain_length),
                holo=bool(row.chain_is_holo),
                uniprot_ids=tuple(
                    str(value) for value in _as_list(row.chain_uniprot_ids)
                ),
            )
            chains[chain.asym_id] = chain
            if chain.is_polypeptide:
                author_to_asym[chain.auth_id] = chain.asym_id
        return chains, author_to_asym

    # Compatibility for in-memory test/custom dataframes that only contain
    # system receptor columns. Such rows can reproduce holo selection but do
    # not contain enough information to identify apo or predicted chains.
    for _, row in entry_rows.iterrows():
        asyms = _as_list(row["system_protein_chains_asym_id"])
        auths = _as_list(row["system_protein_chains_auth_id"])
        lengths = _as_list(row["system_protein_chains_length"])
        for inst_chain, auth, length in zip(asyms, auths, lengths):
            asym = str(inst_chain).split(".", 1)[1]
            if asym not in chains:
                chains[asym] = ChainView(
                    asym_id=asym,
                    auth_id=str(auth),
                    length=int(length),
                )
                author_to_asym[str(auth)] = asym
    return chains, author_to_asym


def _protein_instance_chains(
    instance_chains: Iterable[str], chains: dict[str, ChainView]
) -> list[str]:
    """Filter receptor instance chains to polypeptides known by the chain index."""
    result = []
    for instance_chain in instance_chains:
        value = str(instance_chain)
        asym_id = value.split(".", maxsplit=1)[-1]
        chain = chains.get(asym_id)
        if chain is not None and chain.is_polypeptide:
            result.append(value)
    return sorted(set(result))


def _make_ligand_view(
    row: pd.Series,
    *,
    pdb_id: str,
    system_id: str,
    chains: dict[str, ChainView],
) -> LigandView:
    """Build a ligand view from one published-index row."""
    instance_chain = str(row["ligand_instance_chain"])
    asym_id = str(row["ligand_asym_id"])
    ligand_id = str(row["ligand_id"])
    pocket_n2i: dict[str, dict[int, int]] = defaultdict(dict)
    interactions: dict[str, dict[int, Counter[str]]] = defaultdict(
        lambda: defaultdict(Counter)
    )
    protein_chains = _protein_instance_chains(
        _as_list(row["ligand_protein_chains_asym_id"]), chains
    )
    protein_chain_set = set(protein_chains)
    for value in _as_list(row["ligand_neighboring_residues"]):
        inst, rnum, ridx = _parse_neighboring_residue(value)
        if inst in protein_chain_set:
            pocket_n2i[inst][rnum] = ridx
    for value in _as_list(row["ligand_interactions"]):
        inst, rnum, itype = _parse_interaction(value)
        if inst in protein_chain_set:
            interactions[inst][rnum][itype] += 1
    num_pocket_residues = sum(len(residues) for residues in pocket_n2i.values())
    num_interactions = sum(
        sum(counter.values())
        for residues in interactions.values()
        for counter in residues.values()
    )
    num_unique_interactions = sum(
        len(counter)
        for residues in interactions.values()
        for counter in residues.values()
    )
    if "ligand_is_3d_score_able" in row:
        score_ability = row["ligand_is_3d_score_able"]
        is_3d_score_able = False if pd.isna(score_ability) else bool(score_ability)
    else:
        # V2 indexes predate this annotation; retain their runtime behavior.
        is_3d_score_able = True
    return LigandView(
        id=ligand_id,
        pdb_id=pdb_id,
        system_id=system_id,
        instance_chain=instance_chain,
        asym_id=asym_id,
        is_proper=bool(row["ligand_is_proper"]),
        protein_chains_asym_id=protein_chains,
        num_pocket_residues=num_pocket_residues,
        num_interactions=num_interactions,
        num_unique_interactions=num_unique_interactions,
        is_3d_score_able=is_3d_score_able,
        pocket_residue_number_to_index={k: dict(v) for k, v in pocket_n2i.items()},
        interactions_counter={
            k: {r: Counter(c) for r, c in v.items()} for k, v in interactions.items()
        },
    )


def _make_interface_view(row: pd.Series, *, pdb_id: str) -> InterfaceView:
    """Build a protein-interface view from one annotation row."""
    chain_1_numbers = [
        int(value) for value in _as_list(row["interface_chain_1_residue_numbers"])
    ]
    chain_1_indices = [
        int(value) for value in _as_list(row["interface_chain_1_residue_indices"])
    ]
    chain_2_numbers = [
        int(value) for value in _as_list(row["interface_chain_2_residue_numbers"])
    ]
    chain_2_indices = [
        int(value) for value in _as_list(row["interface_chain_2_residue_indices"])
    ]
    if len(chain_1_numbers) != len(chain_1_indices) or len(chain_2_numbers) != len(
        chain_2_indices
    ):
        raise ValueError(
            f"interface residue mapping lengths differ for {row['system_id']}"
        )
    return InterfaceView(
        id=str(row["system_id"]),
        pdb_id=pdb_id,
        biounit_id=str(row["system_biounit_id"]),
        chain_1=str(row["interface_chain_1"]),
        chain_2=str(row["interface_chain_2"]),
        chain_1_residue_number_to_index=dict(zip(chain_1_numbers, chain_1_indices)),
        chain_2_residue_number_to_index=dict(zip(chain_2_numbers, chain_2_indices)),
        num_contact_residue_pairs=int(row["interface_num_contact_residue_pairs"]),
    )


def entry_views_from_df(
    df: pd.DataFrame,
    *,
    entry_chains: pd.DataFrame | None = None,
    interface_annotations: pd.DataFrame | None = None,
) -> dict[str, EntryView]:
    """Build :class:`EntryView` objects from any DataFrame shaped like the
    published annotation parquet — i.e. one row per
    ``(entry, system, ligand)`` triple. ``interface_annotations`` supplies the
    parallel one-row-per-protein-interface table and permits interface-only
    entries. Pass the one-row-per-chain table to retain receptor
    types and apo/predicted alignment metadata.
    Without it, only protein-only holo chains present on system rows can be
    reconstructed. Annotation-only mixed or nucleic-acid receptors cannot be
    reconstructed because the system-level receptor type cannot be assigned
    to individual chains unambiguously.

    Source-agnostic: works equally on the published parquet read via
    :func:`load_entry_views`, a locally-built parquet, or a freshly
    constructed DataFrame from in-memory ``Entry`` objects
    (``pd.concat([e.to_df() for e in entries.values()])``).
    """
    interface_annotations = (
        interface_annotations
        if interface_annotations is not None
        else pd.DataFrame(columns=["entry_pdb_id"])
    )
    pdb_ids = list(
        dict.fromkeys(
            [str(value) for value in df.get("entry_pdb_id", [])]
            + [str(value) for value in interface_annotations.get("entry_pdb_id", [])]
        )
    )
    views: dict[str, EntryView] = {}
    for pdb_id in pdb_ids:
        entry_rows = df[df["entry_pdb_id"].astype(str) == pdb_id]
        interface_rows = interface_annotations[
            interface_annotations["entry_pdb_id"].astype(str) == pdb_id
        ]
        chain_rows = None
        if entry_chains is not None:
            chain_rows = entry_chains[entry_chains["entry_pdb_id"] == pdb_id]
            if chain_rows.empty:
                receptor_types = {
                    str(value)
                    for value in entry_rows.get(
                        "system_receptor_type", pd.Series(dtype=str)
                    )
                }
                if (
                    not interface_rows.empty
                    or not receptor_types
                    or any(
                        "protein" in receptor_type.split("+")
                        for receptor_type in receptor_types
                    )
                ):
                    raise ValueError(f"No protein chain metadata found for {pdb_id}")
        else:
            receptor_types = {
                str(value)
                for value in entry_rows.get(
                    "system_receptor_type", pd.Series(dtype=str)
                )
                if pd.notna(value)
            }
            unsupported_types = sorted(
                receptor_type
                for receptor_type in receptor_types
                if receptor_type != "protein"
            )
            if unsupported_types:
                raise ValueError(
                    f"entry_chains is required for {pdb_id}: annotation-only "
                    "annotations cannot assign receptor chain types for "
                    f"{unsupported_types}"
                )
        if "system_id" not in entry_rows.columns:
            entry_rows = entry_rows.assign(system_id=pd.Series(dtype="string"))
        chains, author_to_asym = _entry_chains_from_rows(entry_rows, chain_rows)

        systems: dict[str, SystemView] = {}
        for system_id, sys_rows in entry_rows.groupby("system_id", sort=False):
            first = sys_rows.iloc[0]
            ligands: dict[str, LigandView] = {}
            for _, row in sys_rows.iterrows():
                ligand = _make_ligand_view(
                    row,
                    pdb_id=str(pdb_id),
                    system_id=str(system_id),
                    chains=chains,
                )
                ligands[ligand.instance_chain] = ligand
            pocket_n2i: dict[str, dict[int, int]] = defaultdict(dict)
            interactions: dict[str, dict[int, Counter[str]]] = defaultdict(
                lambda: defaultdict(Counter)
            )
            proper_ligands = [ligand for ligand in ligands.values() if ligand.is_proper]
            for ligand in proper_ligands:
                for (
                    instance_chain,
                    pocket_residues,
                ) in ligand.pocket_residue_number_to_index.items():
                    pocket_n2i[instance_chain].update(pocket_residues)
                for (
                    instance_chain,
                    interaction_residues,
                ) in ligand.interactions_counter.items():
                    for residue_number, counter in interaction_residues.items():
                        interactions[instance_chain][residue_number].update(counter)
            systems[system_id] = SystemView(
                id=system_id,
                pdb_id=str(pdb_id),
                system_type=first["system_type"],
                protein_chains_asym_id=_protein_instance_chains(
                    _as_list(first["system_protein_chains_asym_id"]), chains
                ),
                proper_num_pocket_residues=sum(
                    len(residues) for residues in pocket_n2i.values()
                ),
                proper_num_interactions=sum(
                    ligand.num_interactions for ligand in proper_ligands
                ),
                proper_num_unique_interactions=sum(
                    ligand.num_unique_interactions for ligand in proper_ligands
                ),
                receptor_type=str(first.get("system_receptor_type", "protein")),
                pocket_residue_number_to_index={
                    k: dict(v) for k, v in pocket_n2i.items()
                },
                interactions_counter={
                    k: {r: Counter(c) for r, c in v.items()}
                    for k, v in interactions.items()
                },
                ligands=ligands,
            )
        interfaces: dict[str, InterfaceView] = {}
        for _, row in interface_rows.iterrows():
            interface = _make_interface_view(row, pdb_id=pdb_id)
            interfaces[interface.id] = interface
        for interface in interfaces.values():
            missing_chains = {
                chain.split(".", maxsplit=1)[-1]
                for chain in interface.chains
                if chain.split(".", maxsplit=1)[-1] not in chains
            }
            if missing_chains:
                raise ValueError(
                    f"interface {interface.id} references missing chains: "
                    f"{sorted(missing_chains)}"
                )
        views[pdb_id] = EntryView(
            pdb_id=pdb_id,
            chains=chains,
            systems=systems,
            author_to_asym=author_to_asym,
            interfaces=interfaces,
        )
    return views


def load_entry_views(
    *, pdb_ids: Iterable[str], data_dir: Path | None = None
) -> dict[str, EntryView]:
    """Load annotation and chain rows for the requested entries.

    ``data_dir`` selects a local ingest/release root. If omitted, the ligand,
    interface, and chain tables are resolved from the configured PLINDER
    release cache. Production callers should use this loader so all rows come
    from the same release. Direct/custom DataFrames can use
    :func:`entry_views_from_df` instead.
    """
    pdb_ids = sorted(set(pdb_ids))
    if not pdb_ids:
        return {}

    from plinder.core.release import PlinderRelease

    release = PlinderRelease(data_dir)
    if data_dir is None:
        df = query_index(
            columns=["*"],
            filters=[FILTER(("entry_pdb_id", "in", set(pdb_ids)))],
        )
        chain_path = release.fetch("entry_chains")
        interface_path = release.fetch("interface_annotations")
    else:
        annotation_path = release.path("annotation_table")
        if not annotation_path.is_file():
            raise FileNotFoundError(f"missing annotation index: {annotation_path}")
        df = pd.read_parquet(
            annotation_path,
            filters=[("entry_pdb_id", "in", pdb_ids)],
        )
        chain_path = release.path("entry_chains")
        interface_path = release.path("interface_annotations")

    if not chain_path.is_file():
        raise FileNotFoundError(f"missing entry chain index: {chain_path}")
    if not interface_path.is_file():
        raise FileNotFoundError(
            f"missing interface annotation index: {interface_path}"
        )
    entry_chains = pd.read_parquet(
        chain_path,
        filters=[("entry_pdb_id", "in", pdb_ids)],
    )
    interfaces = pd.read_parquet(
        interface_path,
        filters=[("entry_pdb_id", "in", pdb_ids)],
    )
    LOG.info(
        "load_entry_views: %s ligand rows and %s interface rows for %s pdb_ids",
        len(df),
        len(interfaces),
        len(pdb_ids),
    )
    return entry_views_from_df(
        df,
        entry_chains=entry_chains,
        interface_annotations=interfaces,
    )


def load_alignment_entry_views(
    *, lookup_path: Path, pdb_ids: Iterable[str]
) -> dict[str, EntryView]:
    """Load compact chain and selected-residue maps for alignment mapping."""
    selected = sorted(set(pdb_ids))
    if not selected:
        return {}
    frame = pd.read_parquet(
        lookup_path,
        filters=[("entry_pdb_id", "in", selected)],
    )
    views: dict[str, EntryView] = {}
    for pdb_id, rows in frame.groupby("entry_pdb_id", sort=False):
        author_to_asym: dict[str, str] = {}
        pocket_n2i: dict[str, dict[int, int]] = {}
        for row in rows.itertuples(index=False):
            asym_id = str(row.chain_asym_id)
            author_to_asym[str(row.chain_auth_id)] = asym_id
            numbers = [int(value) for value in _as_list(row.selected_residue_numbers)]
            indices = [int(value) for value in _as_list(row.selected_residue_indices)]
            if numbers:
                pocket_n2i[f"1.{asym_id}"] = dict(zip(numbers, indices))
        systems = {}
        if pocket_n2i:
            systems["alignment"] = SystemView(
                id="alignment",
                pdb_id=str(pdb_id),
                system_type="holo",
                protein_chains_asym_id=sorted(pocket_n2i),
                proper_num_pocket_residues=sum(map(len, pocket_n2i.values())),
                proper_num_interactions=0,
                proper_num_unique_interactions=0,
                pocket_residue_number_to_index=pocket_n2i,
            )
        views[str(pdb_id)] = EntryView(
            pdb_id=str(pdb_id),
            chains={},
            systems=systems,
            author_to_asym=author_to_asym,
        )
    LOG.info(
        f"load_alignment_entry_views: {len(frame)} chains for {len(views)} pdb_ids"
    )
    return views
