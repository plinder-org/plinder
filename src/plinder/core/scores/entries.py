# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
"""Lightweight entry/system/chain views loaded from the index parquet.

Consumed by similarity scoring.
"""
from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass, field
from functools import cached_property
from typing import Iterable

from plinder.core.scores.index import query_index
from plinder.core.utils.log import setup_logger

LOG = setup_logger(__name__)


@dataclass(frozen=True)
class ChainView:
    asym_id: str  # label_asym_id (e.g. "A")
    auth_id: str  # auth_asym_id (e.g. "A")
    length: int  # SEQRES length


@dataclass
class SystemView:
    id: str
    pdb_id: str
    system_type: str
    protein_chains_asym_id: list[str]  # instance_chain strings, e.g. ["1.A", "1.B"]
    proper_num_pocket_residues: int
    proper_num_interactions: int
    proper_num_unique_interactions: int
    # instance_chain -> {residue_number: residue_index}
    pocket_residue_number_to_index: dict[str, dict[int, int]] = field(default_factory=dict)
    # instance_chain -> {residue_number: Counter[interaction_type]}
    interactions_counter: dict[str, dict[int, Counter[str]]] = field(default_factory=dict)

    @cached_property
    def pocket_residue_index_to_number(self) -> dict[str, dict[int, int]]:
        return {
            chain: {idx: num for num, idx in mapping.items()}
            for chain, mapping in self.pocket_residue_number_to_index.items()
        }


@dataclass
class EntryView:
    pdb_id: str
    chains: dict[str, ChainView]  # by asym_id
    systems: dict[str, SystemView]  # by system_id
    author_to_asym: dict[str, str]

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

    def chains_for_alignment(self, chain_type: str, aln_type: str) -> list[str]:
        if chain_type != "holo":
            raise NotImplementedError(
                f"chains_for_alignment chain_type={chain_type!r} not supported by "
                "EntryView; apo/pred scoring needs additional columns"
            )
        receptor_asym_ids = {
            instance_chain.split(".")[1]
            for system in self.systems.values()
            if system.system_type == "holo"
            for instance_chain in system.protein_chains_asym_id
        }
        chains = sorted(
            self.chains[asym].auth_id
            for asym in receptor_asym_ids
            if asym in self.chains
        )
        if aln_type == "foldseek":
            return [f"pdb_0000{self.pdb_id}_xyz-enrich_{c}" for c in chains]
        if aln_type == "mmseqs":
            return [f"{self.pdb_id}_{c}" for c in chains]
        raise ValueError(f"unknown aln_type={aln_type!r}")


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


def _as_list(value: object) -> list:
    """Normalize an array-valued cell (numpy array / None / NaN) to a list."""
    if value is None:
        return []
    try:
        return list(value)
    except TypeError:
        return []


def entry_views_from_df(df: "pd.DataFrame") -> dict[str, EntryView]:
    """Build :class:`EntryView` objects from any DataFrame shaped like the
    published index parquet — i.e. one row per ``(entry, system, ligand)``
    triple with the same column names produced by ``Entry.to_df()``.

    Source-agnostic: works equally on the published parquet read via
    :func:`load_entry_views`, a locally-built parquet, or a freshly
    constructed DataFrame from in-memory ``Entry`` objects
    (``pd.concat([e.to_df() for e in entries.values()])``).
    """
    views: dict[str, EntryView] = {}
    for pdb_id, entry_rows in df.groupby("entry_pdb_id", sort=False):
        chains: dict[str, ChainView] = {}
        author_to_asym: dict[str, str] = {}
        for _, row in entry_rows.iterrows():
            asyms = _as_list(row["system_protein_chains_asym_id"])
            auths = _as_list(row["system_protein_chains_auth_id"])
            lengths = _as_list(row["system_protein_chains_length"])
            for inst_chain, auth, length in zip(asyms, auths, lengths):
                asym = inst_chain.split(".", 1)[1]
                if asym not in chains:
                    chains[asym] = ChainView(
                        asym_id=asym, auth_id=auth, length=int(length)
                    )
                    author_to_asym[auth] = asym

        systems: dict[str, SystemView] = {}
        for system_id, sys_rows in entry_rows.groupby("system_id", sort=False):
            first = sys_rows.iloc[0]
            proper_rows = sys_rows[sys_rows["ligand_is_proper"].astype(bool)]
            pocket_n2i: dict[str, dict[int, int]] = defaultdict(dict)
            interactions: dict[str, dict[int, Counter[str]]] = defaultdict(
                lambda: defaultdict(Counter)
            )
            for _, lig in proper_rows.iterrows():
                for s in _as_list(lig["ligand_neighboring_residues"]):
                    inst, rnum, ridx = _parse_neighboring_residue(s)
                    pocket_n2i[inst][rnum] = ridx
                for s in _as_list(lig["ligand_interactions"]):
                    inst, rnum, itype = _parse_interaction(s)
                    interactions[inst][rnum][itype] += 1
            systems[system_id] = SystemView(
                id=system_id,
                pdb_id=str(pdb_id),
                system_type=first["system_type"],
                protein_chains_asym_id=_as_list(
                    first["system_protein_chains_asym_id"]
                ),
                proper_num_pocket_residues=int(first["system_proper_num_pocket_residues"]),
                proper_num_interactions=int(first["system_proper_num_interactions"]),
                proper_num_unique_interactions=int(
                    first["system_proper_num_unique_interactions"]
                ),
                pocket_residue_number_to_index={
                    k: dict(v) for k, v in pocket_n2i.items()
                },
                interactions_counter={
                    k: {r: Counter(c) for r, c in v.items()}
                    for k, v in interactions.items()
                },
            )
        views[str(pdb_id)] = EntryView(
            pdb_id=str(pdb_id),
            chains=chains,
            systems=systems,
            author_to_asym=author_to_asym,
        )
    return views


def load_entry_views(*, pdb_ids: Iterable[str]) -> dict[str, EntryView]:
    """Build :class:`EntryView` objects for the given pdb_ids from the
    published plinder index parquet. Thin convenience wrapper around
    :func:`entry_views_from_df` — for sources other than the published
    index, call ``entry_views_from_df`` directly with your own DataFrame.
    """
    pdb_ids = list(pdb_ids)
    df = query_index(columns=["*"], splits=["*"])
    df = df[df["entry_pdb_id"].isin(pdb_ids)]
    LOG.info(f"load_entry_views: {len(df)} rows for {len(pdb_ids)} pdb_ids")
    return entry_views_from_df(df)
