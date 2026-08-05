# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
from dataclasses import dataclass, field
from functools import partial
from typing import Any, Optional

from omegaconf import DictConfig

from plinder.core.scores.metrics import DEFAULT_CLUSTER_METRICS
from plinder.core.utils import config as _config

METRICS = list(DEFAULT_CLUSTER_METRICS)


@dataclass
class FlowConfig:
    """
    Control re-processing logic for the pipeline

    Note
    ----
    See plinder.core.utils.config.DataConfig for more details

    Attributes
    ----------
    run_specific_stages : str, default=""
        if set, comma-separated list of specific stages to run
    skip_specific_stages : str, default=""
        if set, comma-separated list of specific stages to skip
    download_rcsb_files_batch_size: int
        Target number of two_char_codes per task batch.
    annotation_batch_size : int, default=220
        How many system annotations to generate in a given chunk
    skip_existing_entries : bool, default=True
        if the per-entry annotation parquet already exists, skip generation
    make_entries_cpu : int, default=1
        misguided experiments in multiprocessing over C++ libs (bad idea)
    make_entries_mode : str, default="all"
        Generate ligands and interfaces, only ligands, or only interfaces.
    """

    run_specific_stages: Any = ""
    skip_specific_stages: Any = ""
    # 1060 codes but RCSB limits rsync connections
    download_rcsb_files_batch_size: int = 4

    # 220K PDB IDs but cluster has max workers
    make_entries_batch_size: int = 220
    make_entries_force_update: bool = False
    make_entries_cpu: int = 4
    make_entries_mode: str = "all"

    collate_entries_batch_size: int = 4
    collate_entries_cpu: int = 2
    collate_entries_memory_limit: str = "7GB"
    finalize_entries_cpu: int = 4
    finalize_entries_memory_limit: str = "32GB"

    make_sub_dbs_cpu: int = 4
    make_scorers_cpu: int = 4
    download_alternative_datasets_threads: int = 10
    make_dbs_cpu: int = 4

    make_ligands_batch_size: int = 100
    make_ligands_force_update: bool = False

    cluster_metrics: list[str] = field(default_factory=lambda: METRICS.copy())
    cluster_thresholds: list[int] = field(default_factory=lambda: [30, 50, 70, 90, 100])
    symmetric_edge_source_batch_size: int = 20
    symmetric_edge_bucket_count: int = 64
    component_reduction_source_batch_size: int = 1
    component_reduction_metric_workers: int = 4
    make_communities_cpu: int = 4
    make_components_force_update: bool = True
    make_components_stop_on_cluster: int = 0

    run_batch_searches_batch_size: int = 5_000
    map_batch_alignments_batch_size: int = 25
    make_interface_scores_batch_size: int = 1
    make_interface_scores_cpu: int = 4
    make_interface_scores_memory_limit: str = "32GB"
    make_batch_scores_batch_size: int = 90
    make_batch_scores_cpu: int = 4
    make_batch_scores_force_update: bool = False
    make_ligand_3d_scores_batch_size: int = 30_000
    make_ligand_3d_scores_cpu: int = 1
    collate_ligand_3d_candidates_batch_size: int = 4
    collate_ligand_3d_scores_batch_size: int = 4
    merge_ligand_3d_scores_batch_size: int = 4
    collate_partitions_cpu: int = 4
    collate_partitions_memory_limit: str = "7GB"

    make_links_cpu: int = 8
    make_linked_structures_cpu: int = 8
    make_linked_structures_force_update: bool = False
    score_linked_structures_cpu: int = 8
    score_linked_structures_batch_size: int = 100
    score_linked_structures_force_update: bool = False
    sub_databases: Any = "apo,pred"

    split_config_dir: str = ""

    def __post_init__(self) -> None:
        if self.make_entries_mode not in {"all", "ligands", "interfaces"}:
            raise ValueError(
                "flow.make_entries_mode must be all, ligands, or interfaces"
            )
        if self.collate_entries_batch_size < 1:
            raise ValueError("flow.collate_entries_batch_size must be positive")
        if self.collate_entries_cpu < 1 or self.finalize_entries_cpu < 1:
            raise ValueError("entry collation CPU counts must be positive")
        for name in [
            "component_reduction_metric_workers",
            "symmetric_edge_source_batch_size",
            "symmetric_edge_bucket_count",
            "make_interface_scores_batch_size",
            "make_interface_scores_cpu",
            "make_ligand_3d_scores_batch_size",
            "make_ligand_3d_scores_cpu",
            "collate_ligand_3d_candidates_batch_size",
            "collate_ligand_3d_scores_batch_size",
            "merge_ligand_3d_scores_batch_size",
        ]:
            if getattr(self, name) < 1:
                raise ValueError(f"flow.{name} must be positive")
        if isinstance(self.run_specific_stages, str):
            self.run_specific_stages = [
                stage for stage in self.run_specific_stages.split(",") if stage
            ]
        if isinstance(self.skip_specific_stages, str):
            self.skip_specific_stages = [
                stage for stage in self.skip_specific_stages.split(",") if stage
            ]
        if isinstance(self.sub_databases, str):
            self.sub_databases = [db for db in self.sub_databases.split(",") if db]


@dataclass
class SourceConfig:
    """Locations of the source archives consumed by V3 entry ingest.

    Empty roots use ``PLINDER_PDB_NEXTGEN_ROOT`` and
    ``PLINDER_VALIDATION_ROOT`` when set, then fall back to the Metaflow-local
    ``ingest`` and ``reports`` directories.
    """

    pdb_nextgen_root: str = ""
    validation_root: str = ""
    seqres_path: str = ""
    discovery_threads: int = 8

    def __post_init__(self) -> None:
        if self.discovery_threads < 1:
            raise ValueError("source.discovery_threads must be positive")


@dataclass
class FoldseekConfig:
    alignment_type: int = 2
    score_type: str = "lddt"
    evalue: float = 0.01  # pinder uses default=0.05
    max_seqs: int = 10_000
    sensitivity: float = 11.0  # pinder uses default=11.0
    min_seq_id: float = 0.2
    coverage: float = 0.0
    alignment_filename: str = "alignment.txt"

    def __post_init__(self) -> None:
        if self.max_seqs < 1:
            raise ValueError("foldseek.max_seqs must be positive")
        if not 0 <= self.min_seq_id <= 1:
            raise ValueError("foldseek.min_seq_id must be in [0, 1]")
        for attr, allowed in [
            ("alignment_type", [1, 2]),
            ("score_type", ["lddt", "alntmscore"]),
        ]:
            passed = getattr(self, attr, None)
            if passed not in allowed:  # type: ignore
                raise ValueError(
                    f"{self.__class__.__name__}.{attr} must be in {allowed}"
                )


@dataclass
class MMSeqsConfig:
    score_type: str = "pident"
    evalue: float = 0.01  # pinder uses default=0.05
    max_seqs: int = 10_000
    sensitivity: float = 11.0  # pinder uses default=11.0
    min_seq_id: float = 0.2
    coverage: float = 0.0
    alignment_filename: str = "alignment.txt"

    def __post_init__(self) -> None:
        if self.max_seqs < 1:
            raise ValueError("mmseqs.max_seqs must be positive")
        if not 0 <= self.min_seq_id <= 1:
            raise ValueError("mmseqs.min_seq_id must be in [0, 1]")
        for attr, allowed in [
            ("score_type", ["pident"]),
        ]:
            if getattr(self, attr) not in allowed:
                raise ValueError(
                    f"{self.__class__.__name__}.{attr} must be in {allowed}"
                )


@dataclass
class GraphConfig:
    pass


@dataclass
class ScorerConfig:
    wipe_partition: bool = False
    rerun_existing_batch: bool = False
    minimum_threshold: float = 0.3
    minimum_thresholds: dict[str, float] = field(default_factory=dict)
    max_alignment_rows_per_query: int = 5_000_000
    max_query_protein_chains: int = 30
    max_query_proper_ligand_chains: int = 30
    sub_databases: Any = "holo,apo,pred"

    def __post_init__(self) -> None:
        if isinstance(self.sub_databases, str):
            self.sub_databases = [db for db in self.sub_databases.split(",") if db]
        allowed_dbs = ["holo", "apo", "pred"]
        for db in self.sub_databases:
            if db not in allowed_dbs:
                raise ValueError(
                    f"{self.__class__.__name__}.sub_databases must be in {allowed_dbs}"
                )
        if self.minimum_threshold < 0 or self.minimum_threshold > 1:
            raise ValueError(
                f"{self.__class__.__name__}.minimum_threshold must be in [0, 1]"
            )
        for metric, threshold in self.minimum_thresholds.items():
            if threshold < 0 or threshold > 1:
                raise ValueError(
                    f"scorer.minimum_thresholds[{metric!r}] must be in [0, 1]"
                )
        if self.max_alignment_rows_per_query < 1:
            raise ValueError("scorer.max_alignment_rows_per_query must be positive")
        if self.max_query_protein_chains < 1 or self.max_query_proper_ligand_chains < 1:
            raise ValueError("scorer query system chain limits must be positive")


@dataclass
class EntryConfig:
    # TODO-tjd: deduplicate with AnnotationConfig
    neighboring_residue_threshold: float = 6.0
    neighboring_ligand_threshold: float = 4.0
    min_polymer_size: int = 12
    plip_complex_threshold: float = 10.0
    min_shared_pocket_members: int = 3
    data_dir: Optional[str] = None
    save_folder: Optional[str] = None


@dataclass
class AnnotationConfig:
    """
    Configuration model for tuning how Entry.from_cif_file
    generates systems and entry metadata.

    Note
    ----
    See plinder.data.get_system_annotations.GetPlinderAnnotation
    for more details. This purposely omits the cif_file Path argument
    because that will necessarily be provided by consuming code.
    """

    neighboring_residue_threshold: float = 6.0
    neighboring_ligand_threshold: float = 4.0
    min_polymer_size: int = 12
    min_shared_pocket_members: int = 3


@dataclass
class InterfaceConfig:
    """Protein-interface definitions used during entry ingest."""

    contact_radius: float = 10.0
    min_chain_length: int = 12
    min_interface_residues: int = 7
    annotate_prodigy: bool = True

    def __post_init__(self) -> None:
        if self.contact_radius <= 0:
            raise ValueError("interface.contact_radius must be positive")
        if self.min_chain_length < 1:
            raise ValueError("interface.min_chain_length must be positive")
        if self.min_interface_residues < 1:
            raise ValueError("interface.min_interface_residues must be positive")


""" From
OleinikovasV
OleinikovasV commented Apr 22, 2024

Added updated artifacts list and curation, please, review the logic!

"""


@dataclass
class LigandConfig:
    minimum_similarity: float = 30.0
    cofactor_similarity_threshold: float = 90.0
    number_id_col: str = "ligand_smiles_id"

    def __post_init__(self) -> None:
        for name in [
            "minimum_similarity",
            "cofactor_similarity_threshold",
        ]:
            value = float(getattr(self, name))
            if not 0 <= value <= 100:
                raise ValueError(f"ligand.{name} must be between 0 and 100")
        if self.minimum_similarity > 90:
            raise ValueError(
                "ligand.minimum_similarity must not exceed the fixed 90-percent "
                "frequency cluster threshold"
            )


SCHEMA = {
    "flow": FlowConfig,
    "source": SourceConfig,
    "foldseek": FoldseekConfig,
    "mmseqs": MMSeqsConfig,
    "graph": GraphConfig,
    "annotation": AnnotationConfig,
    "interface": InterfaceConfig,
    "entry": EntryConfig,
    "scorer": ScorerConfig,
    "ligand": LigandConfig,
}

SCHEMA.update(_config.SCHEMA)

_get_config = partial(_config._config, schema=SCHEMA, package_schema="data")


def get_config(**kwargs: Any) -> DictConfig:
    """Load and cross-validate the data-pipeline configuration."""
    cfg = _get_config(**kwargs)
    unsupported_metrics = sorted(
        set(cfg.flow.cluster_metrics).difference(DEFAULT_CLUSTER_METRICS)
    )
    if unsupported_metrics:
        raise ValueError(
            "unsupported release clustering metrics: "
            f"{unsupported_metrics}; use sucos_shape_pocket_qcov instead of "
            "raw shape, color, or SuCOS"
        )
    required_thresholds = [90.0]
    if "tanimoto_similarity_ecfp4_1024" in cfg.flow.cluster_metrics:
        required_thresholds.extend(
            float(value) for value in cfg.flow.cluster_thresholds
        )
    lowest_required_threshold = min(required_thresholds)
    if cfg.ligand.minimum_similarity > lowest_required_threshold:
        raise ValueError(
            "ligand.minimum_similarity must not exceed the lowest requested "
            "Tanimoto clustering threshold "
            f"({lowest_required_threshold:g})"
        )
    derived_metrics = set(cfg.flow.cluster_metrics).difference(
        {"tanimoto_similarity_ecfp4_1024"}
    )
    lowest_cluster_threshold = min(
        float(value) for value in cfg.flow.cluster_thresholds
    )
    if (
        derived_metrics
        and float(cfg.scorer.minimum_threshold) * 100 > lowest_cluster_threshold
    ):
        raise ValueError(
            "scorer.minimum_threshold must not exceed the lowest requested "
            "derived-score clustering threshold "
            f"({lowest_cluster_threshold:g})"
        )
    return cfg
