# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
"""Paths and table keys for a published PLINDER release."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

from plinder.core.utils import cpl
from plinder.core.utils.config import get_config

RELEASE_PATHS = {
    "annotation_table": "index/annotation_table.parquet",
    "entry_chains": "index/entry_chains.parquet",
    "entry_biounit_chains": "index/entry_biounit_chains.parquet",
    "entry_metadata": "index/entry_metadata.parquet",
    "entry_sources": "index/entry_sources.parquet",
    "linked_apo_structures": "index/linked_apo_structures.parquet",
    "interface_annotations": "index/interface_annotation_table.parquet",
    "alignment_chain_lookup": "index/alignment_chain_lookup.parquet",
    "ligand_pocket_membership": "index/ligand_pocket_membership.parquet",
    "ligand_pocket_representatives": ("index/ligand_pocket_representatives.parquet"),
    "interface_half_representatives": ("index/interface_half_representatives.parquet"),
    "interface_membership": "index/interface_membership.parquet",
    "interface_representatives": "index/interface_representatives.parquet",
    "alignments": "alignments",
    "alignment_shard": (
        "alignments/search_db={search_db}/alignment_type={alignment_type}/"
        "shard={shard}.parquet"
    ),
    "ligand_archives": "ligand_archives",
    "ligand_archive": "ligand_archives/{shard}.parquet",
    "ligand_archives_manifest": "ligand_archives/manifest.json",
    "ligand_scores": "ligand_scores",
    "interface_scores": "interface_scores",
    "ligand_clusters": "ligand_clusters",
    "ligand_sampling": "ligand_sampling",
    "interface_clusters": "interface_clusters",
    "interface_sampling": "interface_sampling",
    "search_databases": "search_databases",
    "search_database": "search_databases/holo_{backend}",
    "sucos_export": "exports/all_sucos_shape_pocket_qcov.parquet",
    "exports": "exports",
}


RELEASE_TABLES: dict[str, dict[str, Any]] = {
    "annotation": {
        "artifact": "annotation_table",
        "row_grain": "ligand",
        "primary_key": ("ligand_id",),
    },
    "entry_chains": {
        "artifact": "entry_chains",
        "row_grain": "polymer chain in a PDB entry",
        "primary_key": ("entry_pdb_id", "chain_asym_id"),
    },
    "entry_biounit_chains": {
        "artifact": "entry_biounit_chains",
        "row_grain": "chain instance in a biological assembly",
        "primary_key": ("entry_pdb_id", "biounit_id", "chain_instance"),
    },
    "entry_metadata": {
        "artifact": "entry_metadata",
        "row_grain": "PDB entry",
        "primary_key": ("entry_pdb_id",),
    },
    "entry_sources": {
        "artifact": "entry_sources",
        "row_grain": "PDB entry",
        "primary_key": ("entry_pdb_id",),
    },
    "interface_annotations": {
        "artifact": "interface_annotations",
        "row_grain": "protein-chain interface in a biological assembly",
        "primary_key": ("system_id",),
    },
    "alignment_chain_lookup": {
        "artifact": "alignment_chain_lookup",
        "row_grain": "protein chain in the release search databases",
        "primary_key": ("entry_pdb_id", "chain_asym_id"),
    },
    "ligand_pocket_membership": {
        "artifact": "ligand_pocket_membership",
        "row_grain": "ligand directed-set-cover assignment",
        "primary_key": ("ligand_id",),
    },
    "ligand_pocket_representatives": {
        "artifact": "ligand_pocket_representatives",
        "row_grain": "ligand-pocket representative",
        "primary_key": ("representative_ligand_id",),
    },
    "interface_half_representatives": {
        "artifact": "interface_half_representatives",
        "row_grain": "representative half-interface",
        "primary_key": ("half_interface_id",),
    },
    "interface_membership": {
        "artifact": "interface_membership",
        "row_grain": "protein-interface representative assignment",
        "primary_key": ("system_id",),
    },
    "interface_representatives": {
        "artifact": "interface_representatives",
        "row_grain": "protein-interface representative",
        "primary_key": ("representative_system_id",),
    },
}


@dataclass(frozen=True)
class PlinderRelease:
    """Resolve named paths within a local or configured PLINDER release."""

    data_dir: Path | None = None

    def _relative_path(self, name: str, **parameters: str) -> Path:
        try:
            template = RELEASE_PATHS[name]
        except KeyError as exc:
            choices = ", ".join(sorted(RELEASE_PATHS))
            raise KeyError(
                f"unknown release artifact {name!r}; choose from: {choices}"
            ) from exc

        for parameter, value in parameters.items():
            value = str(value)
            if not value or value in {".", ".."} or Path(value).name != value:
                raise ValueError(
                    f"artifact parameter {parameter!r} must be one path component"
                )
        try:
            return Path(template.format(**parameters))
        except KeyError as exc:
            raise ValueError(f"missing artifact parameter {exc.args[0]!r}") from exc

    def path(self, name: str, **parameters: str) -> Path:
        """Return the expected local path without fetching it."""
        root = self.data_dir or Path(get_config().data.plinder_dir)
        return Path(root) / self._relative_path(name, **parameters)

    def fetch(self, name: str, **parameters: str) -> Path:
        """Fetch a configured artifact or check an explicit local release."""
        relative = self._relative_path(name, **parameters)
        path = (
            Path(self.data_dir) / relative
            if self.data_dir is not None
            else cpl.get_plinder_path(rel=relative.as_posix())
        )
        if not path.exists():
            raise FileNotFoundError(
                f"PLINDER release artifact {name!r} is unavailable at {path}"
            )
        return path

    def table(self, name: str) -> dict[str, Any]:
        """Return the path name, row grain, and key for a release table."""
        try:
            return RELEASE_TABLES[name]
        except KeyError as exc:
            choices = ", ".join(sorted(RELEASE_TABLES))
            raise KeyError(
                f"unknown release table {name!r}; choose from: {choices}"
            ) from exc
