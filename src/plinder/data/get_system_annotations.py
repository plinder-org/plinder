# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
from __future__ import annotations

import logging
from pathlib import Path
from typing import Any, Dict, Optional

import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

from plinder.core.utils.log import setup_logger
from plinder.data.annotations.aggregate_annotations import Entry
from plinder.data.annotations.interface_utils import (
    DEFAULT_MIN_INTERFACE_RESIDUES,
    protein_interfaces_to_table,
)

LOG = setup_logger(__name__, log_level=logging.DEBUG)


class GetPlinderAnnotation:
    def __init__(
        self,
        mmcif_file: Path,
        validation_xml: Path,
        save_folder: Optional[Path] = None,
        data_dir: Optional[Path] = None,
        neighboring_residue_threshold: float = 6.0,
        neighboring_ligand_threshold: float = 4.0,
        min_polymer_size: int = 12,
        min_shared_pocket_members: int = 3,
        symmetry_mate_contact_threshold: float = 5.0,
        entry_cfg: Optional[Dict[Any, Any]] = None,
        interface_cfg: Optional[Dict[Any, Any]] = None,
    ) -> None:
        self.mmcif_file = mmcif_file
        self.validation_xml = Path(validation_xml)
        self.save_folder = Path(save_folder) if save_folder is not None else None
        self.data_dir = Path(data_dir) if data_dir is not None else None
        self.neighboring_residue_threshold = neighboring_residue_threshold
        self.neighboring_ligand_threshold = neighboring_ligand_threshold
        self.min_polymer_size = min_polymer_size
        self.min_shared_pocket_members = min_shared_pocket_members
        self.symmetry_mate_contact_threshold = symmetry_mate_contact_threshold
        self.entry_cfg = entry_cfg
        self.interface_cfg = interface_cfg

    def _entry_options(
        self, *, include_ligands: bool, include_interfaces: bool
    ) -> dict[str, Any]:
        """Resolve shared entry and protein-interface options."""
        entry_cfg: dict[str, Any] = dict(
            neighboring_residue_threshold=self.neighboring_residue_threshold,
            neighboring_ligand_threshold=self.neighboring_ligand_threshold,
            min_polymer_size=self.min_polymer_size,
            min_shared_pocket_members=self.min_shared_pocket_members,
            save_folder=self.save_folder if include_ligands else None,
            data_dir=self.data_dir if include_ligands else None,
            symmetry_mate_contact_threshold=self.symmetry_mate_contact_threshold,
            include_ligands=include_ligands,
            include_interfaces=include_interfaces,
        )
        if self.entry_cfg is not None:
            entry_cfg.update(self.entry_cfg)
        entry_cfg["include_ligands"] = include_ligands
        entry_cfg["include_interfaces"] = include_interfaces
        entry_cfg["save_folder"] = self.save_folder if include_ligands else None
        entry_cfg["data_dir"] = self.data_dir if include_ligands else None
        interface_cfg = dict(self.interface_cfg or {})
        if interface_cfg:
            entry_cfg.update(
                {
                    "interface_contact_radius": interface_cfg.get(
                        "contact_radius", 10.0
                    ),
                    "interface_min_chain_length": interface_cfg.get(
                        "min_chain_length", 12
                    ),
                    "interface_min_residues": interface_cfg.get(
                        "min_interface_residues", DEFAULT_MIN_INTERFACE_RESIDUES
                    ),
                    "interface_annotate_prodigy": interface_cfg.get(
                        "annotate_prodigy", True
                    ),
                }
            )
        return entry_cfg

    def _interface_table(self, entry_cfg: dict[str, Any]) -> pa.Table:
        return protein_interfaces_to_table(
            self.entry.interfaces,
            min_interface_residues=int(
                entry_cfg.get("interface_min_residues", DEFAULT_MIN_INTERFACE_RESIDUES)
            ),
        )

    def _write_shared_sidecars(
        self,
        entry_folder: Path,
        interface_table: pa.Table,
        *,
        replace_interfaces: bool = True,
        preserve_existing_shared: bool = False,
    ) -> None:
        """Write a new entry's sidecars or extend the preserved shared tables."""
        from plinder.data.annotations.cif_utils import (
            get_mmcif_revision,
            read_mmcif_container,
        )

        entry_folder.mkdir(parents=True, exist_ok=True)

        def write_dataframe(
            path: Path,
            frame: pd.DataFrame,
            *,
            merge_keys: tuple[str, ...] | None = None,
        ) -> None:
            if preserve_existing_shared and path.is_file():
                if merge_keys is None:
                    return
                existing = pd.read_parquet(path)
                if not set(merge_keys).issubset(existing.columns):
                    missing = sorted(set(merge_keys).difference(existing.columns))
                    raise ValueError(f"{path} is missing merge keys: {missing}")
                if not set(merge_keys).issubset(frame.columns):
                    missing = sorted(set(merge_keys).difference(frame.columns))
                    raise ValueError(
                        f"incoming {path} is missing merge keys: {missing}"
                    )
                new_columns = [
                    column for column in frame.columns if column not in existing.columns
                ]
                if new_columns:
                    additions = frame.loc[:, [*merge_keys, *new_columns]]
                    existing = existing.merge(
                        additions,
                        on=list(merge_keys),
                        how="left",
                        validate="one_to_one",
                    )
                # Preserve ligand-derived values for existing rows while
                # adding newly available columns and protein chains discovered
                # by interface-only ingest.
                incoming = frame.reindex(columns=existing.columns)
                existing_keys = pd.MultiIndex.from_frame(
                    existing.loc[:, list(merge_keys)]
                )
                incoming_keys = pd.MultiIndex.from_frame(
                    incoming.loc[:, list(merge_keys)]
                )
                missing_rows = incoming.loc[~incoming_keys.isin(existing_keys)]
                if missing_rows.empty and not new_columns:
                    return
                frame = pd.concat([existing, missing_rows], ignore_index=True)
            temporary = path.with_suffix(".parquet.tmp")
            frame.to_parquet(temporary, index=False)
            temporary.replace(path)

        write_dataframe(
            entry_folder / "entry_chains.parquet",
            self.entry.chains_to_df(),
            merge_keys=("entry_pdb_id", "chain_asym_id"),
        )
        write_dataframe(
            entry_folder / "entry_metadata.parquet",
            self.entry.metadata_to_df(),
            merge_keys=("entry_pdb_id",),
        )
        interface_path = entry_folder / "interfaces.parquet"
        if replace_interfaces or not interface_path.is_file():
            temporary = interface_path.with_suffix(".parquet.tmp")
            pq.write_table(interface_table, temporary)
            temporary.replace(interface_path)
        write_dataframe(
            entry_folder / "entry_biounit_chains.parquet",
            self.entry.biounit_chains_to_df(),
            merge_keys=("entry_pdb_id", "biounit_id", "chain_instance"),
        )
        source_path = entry_folder / "entry_source.parquet"
        if not preserve_existing_shared or not source_path.is_file():
            major_revision, minor_revision = get_mmcif_revision(
                read_mmcif_container(self.mmcif_file)
            )
            write_dataframe(
                source_path,
                pd.DataFrame(
                    {
                        "entry_pdb_id": [self.entry.pdb_id],
                        "source_mmcif_major_revision": [major_revision],
                        "source_mmcif_minor_revision": [minor_revision],
                    }
                ),
            )

    def annotate(self, *, include_interfaces: bool = True) -> Optional[pd.DataFrame]:
        """Annotate ligand systems, optionally including protein interfaces."""
        entry_cfg = self._entry_options(
            include_ligands=True,
            include_interfaces=include_interfaces,
        )
        self.entry = Entry.from_cif_file(
            self.mmcif_file,
            **entry_cfg,
        )
        LOG.info(f"created entry for {self.mmcif_file}")
        self.entry.set_validation(self.validation_xml, self.mmcif_file)
        interface_table = self._interface_table(entry_cfg)
        self.interface_df = interface_table.to_pandas()
        self.entry_metadata_df = self.entry.metadata_to_df()
        resolved_save_folder = entry_cfg.get("save_folder")
        if resolved_save_folder is not None:
            if not isinstance(resolved_save_folder, (str, Path)):
                raise TypeError("entry save_folder must be a path")
            entry_folder = Path(resolved_save_folder) / self.entry.pdb_id
            self._write_shared_sidecars(
                entry_folder,
                interface_table,
                replace_interfaces=include_interfaces,
            )
        if not self.entry.systems and not self.entry.interfaces:
            LOG.info(f"no ligand or interface systems for {self.mmcif_file}")
            return None
        self.annotated_df = self.entry.to_df()
        return self.annotated_df

    def annotate_interfaces(self) -> pa.Table:
        """Annotate protein interfaces without deriving or replacing ligands."""
        entry_cfg = self._entry_options(
            include_ligands=False,
            include_interfaces=True,
        )
        self.entry = Entry.from_cif_file(self.mmcif_file, **entry_cfg)
        LOG.info(f"created interface-only entry for {self.mmcif_file}")
        interface_table = self._interface_table(entry_cfg)
        self.interface_df = interface_table.to_pandas()
        if self.save_folder is None:
            return interface_table

        entry_folder = Path(self.save_folder) / self.entry.pdb_id
        ligand_annotation = Path(self.save_folder) / f"{self.entry.pdb_id}.parquet"
        if entry_folder.exists():
            self._write_shared_sidecars(
                entry_folder,
                interface_table,
                preserve_existing_shared=True,
            )
        elif ligand_annotation.is_file():
            raise FileNotFoundError(
                f"ligand annotation exists without its entry sidecars: {ligand_annotation}"
            )
        elif interface_table.num_rows:
            self._write_shared_sidecars(entry_folder, interface_table)
        return interface_table
