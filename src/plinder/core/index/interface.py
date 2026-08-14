# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
"""Public access to one protein interface in a PLINDER release."""

from __future__ import annotations

from functools import cached_property
from pathlib import Path
from typing import TYPE_CHECKING

import biotite.structure as struc
import numpy as np
import pandas as pd

from plinder.core.index.query import query_table
from plinder.core.release import PlinderRelease
from plinder.core.utils.config import get_config
from plinder.core.utils.io import get_pdb_mmcif
from plinder.data.annotations.save_utils import (
    reconstruct_interface,
    save_reconstructed_interface,
)

if TYPE_CHECKING:
    from numpy.typing import NDArray


class PlinderInterface:
    """One annotated protein-chain interface in a biological assembly."""

    def __init__(
        self,
        *,
        system_id: str,
        release: PlinderRelease | None = None,
        source_mmcif: Path | str | None = None,
        reconstruction_dir: Path | str | None = None,
    ) -> None:
        self.system_id = system_id
        self.release = release or PlinderRelease()
        self.source_mmcif = Path(source_mmcif) if source_mmcif is not None else None
        self.reconstruction_dir = (
            Path(reconstruction_dir)
            if reconstruction_dir is not None
            else Path(get_config().data.plinder_dir)
            / "reconstructed_interfaces"
            / system_id
        )
        self._annotation: pd.Series | None = None

    def __repr__(self) -> str:
        return f"PlinderInterface(system_id={self.system_id!r})"

    @property
    def annotation(self) -> pd.Series:
        """Return the unique published annotation row for this interface."""
        if self._annotation is None:
            rows = query_table(
                "interface_annotations",
                filters=[("system_id", "==", self.system_id)],
                release=self.release,
            )
            if rows.empty:
                raise ValueError(f"interface {self.system_id!r} is not in the release")
            if len(rows) != 1:
                raise ValueError(
                    f"interface {self.system_id!r} has {len(rows)} annotation rows"
                )
            self._annotation = rows.iloc[0]
        return self._annotation

    @property
    def entry_pdb_id(self) -> str:
        """Return the source PDB entry ID."""
        return str(self.annotation["entry_pdb_id"])

    @property
    def biounit_id(self) -> str:
        """Return the deposited biological-assembly ID."""
        return str(self.annotation["system_biounit_id"])

    @property
    def chains(self) -> tuple[str, str]:
        """Return the two canonical assembly-chain instance IDs."""
        return (
            str(self.annotation["interface_chain_1"]),
            str(self.annotation["interface_chain_2"]),
        )

    @property
    def source_mmcif_path(self) -> Path:
        """Return the explicit or release-cached deposited PDB mmCIF."""
        if self.source_mmcif is None:
            manifest = self.release.fetch("entry_sources")
            cache_dir = (
                Path(self.release.data_dir) / get_config().data.source_mmcifs
                if self.release.data_dir is not None
                else None
            )
            self.source_mmcif = get_pdb_mmcif(
                self.entry_pdb_id,
                cache_dir=cache_dir,
                manifest_path=manifest,
            )
        if not self.source_mmcif.is_file():
            raise FileNotFoundError(self.source_mmcif)
        return self.source_mmcif

    @cached_property
    def atom_array(self) -> struc.AtomArray:
        """Return the exact biological-assembly chain pair in memory."""
        return reconstruct_interface(self.source_mmcif_path, self.annotation)

    @cached_property
    def chain_structures(self) -> dict[str, struc.AtomArray]:
        """Return the two interface sides keyed by assembly-chain instance."""
        return {
            chain_id: self.atom_array[self.atom_array.chain_id.astype(str) == chain_id]
            for chain_id in self.chains
        }

    @cached_property
    def sequences(self) -> dict[str, str]:
        """Return full polymer sequences keyed by assembly-chain instance."""
        asym_ids = {chain_id.split(".", maxsplit=1)[-1] for chain_id in self.chains}
        rows = query_table(
            "entry_chains",
            columns=["chain_asym_id", "chain_sequence"],
            filters=[
                ("entry_pdb_id", "==", self.entry_pdb_id),
                ("chain_asym_id", "in", sorted(asym_ids)),
            ],
            release=self.release,
        )
        asym_to_sequence = dict(
            zip(rows["chain_asym_id"].astype(str), rows["chain_sequence"].astype(str))
        )
        missing = sorted(asym_ids.difference(asym_to_sequence))
        if missing:
            raise ValueError(
                f"interface {self.system_id!r} has no sequences for chains {missing}"
            )
        return {
            chain_id: asym_to_sequence[chain_id.split(".", maxsplit=1)[-1]]
            for chain_id in self.chains
        }

    @cached_property
    def interface_residue_masks(self) -> dict[str, "NDArray[np.bool_]"]:
        """Return atom masks for the annotated interface residues on each side."""
        masks: dict[str, NDArray[np.bool_]] = {}
        for side, chain_id in enumerate(self.chains, start=1):
            residue_indices = np.asarray(
                self.annotation[f"interface_chain_{side}_residue_indices"],
                dtype=int,
            )
            chain_atom_indices = np.flatnonzero(
                self.atom_array.chain_id.astype(str) == chain_id
            )
            chain = self.atom_array[chain_atom_indices]
            starts = struc.get_residue_starts(chain, add_exclusive_stop=True)
            residue_count = len(starts) - 1
            if np.any(residue_indices < 0) or np.any(residue_indices >= residue_count):
                raise ValueError(
                    f"interface residue indices for {chain_id} fall outside its "
                    f"{residue_count} resolved residues"
                )
            mask = np.zeros(self.atom_array.array_length(), dtype=bool)
            for residue_index in residue_indices:
                start = starts[residue_index]
                stop = starts[residue_index + 1]
                mask[chain_atom_indices[start:stop]] = True
            masks[chain_id] = mask
        return masks

    @cached_property
    def interface_structure(self) -> struc.AtomArray:
        """Return atoms belonging to annotated interface residues."""
        mask = np.zeros(self.atom_array.array_length(), dtype=bool)
        for side_mask in self.interface_residue_masks.values():
            mask |= side_mask
        return self.atom_array[mask]

    def reconstruct(
        self,
        *,
        output_cif: Path | str | None = None,
        overwrite: bool = False,
    ) -> Path:
        """Write the two-chain biological-assembly interface as an mmCIF."""
        if output_cif is None:
            output_cif = self.reconstruction_dir / "interface.cif"
        return save_reconstructed_interface(
            self.source_mmcif_path,
            self.annotation,
            output_cif=output_cif,
            overwrite=overwrite,
            reconstructed=self.atom_array,
        )

    @property
    def interface_cif(self) -> Path:
        """Return the standard reconstructed interface mmCIF path."""
        path = self.reconstruction_dir / "interface.cif"
        if not path.is_file():
            self.reconstruct(output_cif=path)
        return path


__all__ = ["PlinderInterface"]
