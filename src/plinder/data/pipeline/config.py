# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
from dataclasses import dataclass, field
from functools import partial
from typing import Any, Optional

from plinder.core.utils import config as _config

METRICS = [
    # pli_qcov:
    "pli_qcov",
    "pli_unique_qcov",
    # seq_sim:
    "protein_seqsim_qcov_max",
    "protein_seqsim_qcov_weighted_max",
    "protein_seqsim_qcov_weighted_sum",
    "protein_seqsim_max",
    "protein_seqsim_weighted_max",
    "protein_seqsim_weighted_sum",
    # protein
    "protein_fident_qcov_max",
    "protein_fident_qcov_weighted_max",
    "protein_fident_qcov_weighted_sum",
    "protein_fident_max",
    "protein_fident_weighted_max",
    "protein_fident_weighted_sum",
    "protein_lddt_max",
    "protein_lddt_qcov_max",
    "protein_lddt_qcov_weighted_max",
    "protein_lddt_qcov_weighted_sum",
    "protein_lddt_weighted_max",
    "protein_lddt_weighted_sum",
    "protein_qcov_max",
    "protein_qcov_weighted_max",
    "protein_qcov_weighted_sum",
    # pocket
    "pocket_fident_qcov",
    "pocket_fident",
    "pocket_lddt_qcov",
    "pocket_lddt",
    "pocket_qcov",
]


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
        if the entry JSON already exists, skip generation
    make_entries_cpu : int, default=1
        misguided experiments in multiprocessing over C++ libs (bad idea)
    """

    run_specific_stages: Any = ""
    skip_specific_stages: Any = ""
    # 1060 codes but RCSB limits rsync connections
    download_rcsb_files_batch_size: int = 4

    # 220K PDB IDs but cluster has max workers
    make_entries_batch_size: int = 220
    make_entries_force_update: bool = False
    make_entries_cpu: int = 4

    structure_qc_force_update: bool = False

    make_sub_dbs_cpu: int = 4
    make_scorers_cpu: int = 4
    download_alternative_datasets_threads: int = 10

    make_ligands_batch_size: int = 100
    make_ligands_force_update: bool = False

    cluster_metrics: list[str] = field(default_factory=lambda: METRICS.copy())
    cluster_thresholds: list[int] = field(default_factory=lambda: [50, 70, 95, 100])
    make_components_force_update: bool = True
    make_components_stop_on_cluster: int = 0

    run_batch_searches_batch_size: int = 10050
    make_batch_scores_batch_size: int = 90
    make_batch_scores_force_update: bool = False

    make_links_cpu: int = 8
    make_linked_structures_cpu: int = 8
    make_linked_structures_force_update: bool = False
    score_linked_structures_cpu: int = 8
    score_linked_structures_batch_size: int = 100
    score_linked_structures_force_update: bool = False
    sub_databases: Any = "apo,pred"

    split_config_dir: str = ""
    test_leakage: bool = False

    def __post_init__(self) -> None:
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
class FoldseekConfig:
    alignment_type: int = 2
    score_type: str = "lddt"
    evalue: float = 0.01  # pinder uses default=0.05
    max_seqs: int = 5000  # pinder uses default=1000
    sensitivity: float = 11.0  # pinder uses default=11.0
    min_seq_id: float = 0.2
    coverage: float = 0.0
    alignment_filename: str = "alignment.txt"

    def __post_init__(self) -> None:
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
    max_seqs: int = 5000  # pinder uses default=1000
    sensitivity: float = 11.0  # pinder uses default=11.0
    min_seq_id: float = 0.2
    coverage: float = 0.0
    alignment_filename: str = "alignment.txt"

    def __post_init__(self) -> None:
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
    minimum_threshold: float = 0.2
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


@dataclass
class EntryConfig:
    # TODO-tjd: deduplicate with AnnotationConfig
    max_protein_chains_to_save: int = 5
    max_ligand_chains_to_save: int = 5
    neighboring_residue_threshold: float = 6.0
    neighboring_ligand_threshold: float = 4.0
    min_polymer_size: int = 10
    max_non_small_mol_ligand_length: int = 20
    plip_complex_threshold: float = 10.0
    save_folder: Optional[str] = None
    skip_save_systems: bool = False


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
    min_polymer_size: int = 10


""" From
OleinikovasV
OleinikovasV commented Apr 22, 2024

Added updated artifacts list and curation, please, review the logic!

"""


@dataclass
class LigandConfig:
    radius: int = 2
    nbits: int = 1024
    ligand_id_split_char: str = "__"
    save_top_k_similar_ligands: int = 5000
    multiply_by: int = 100
    number_id_col: str = "number_id_by_inchikeys"
    score_name: str = "tanimoto_similarity_max"


SCHEMA = {
    "flow": FlowConfig,
    "foldseek": FoldseekConfig,
    "mmseqs": MMSeqsConfig,
    "graph": GraphConfig,
    "annotation": AnnotationConfig,
    "entry": EntryConfig,
    "scorer": ScorerConfig,
    "ligand": LigandConfig,
}

SCHEMA.update(_config.SCHEMA)

get_config = partial(_config._config, schema=SCHEMA, package_schema="data")
