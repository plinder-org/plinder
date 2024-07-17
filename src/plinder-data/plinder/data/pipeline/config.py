# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
import hashlib
import json
from dataclasses import dataclass, field
from functools import partial
from io import StringIO
from os import getenv
from pathlib import Path
from typing import Any, Optional
from dataclasses import dataclass, field

from omegaconf import DictConfig, ListConfig, OmegaConf

from plinder.core.utils import config as _config


METRICS = [
    # pli_qcov:
    "pli_qcov_max",
    "pli_qcov_weighted_max",
    "pli_qcov",
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
    "pocket_fident_qcov_max",
    "pocket_fident_qcov_weighted_max",
    "pocket_fident_qcov",
    "pocket_fident_max",
    "pocket_fident_weighted_max",
    "pocket_fident",
    "pocket_lddt_max",
    "pocket_lddt_qcov_max",
    "pocket_lddt_qcov_weighted_max",
    "pocket_lddt_qcov",
    "pocket_lddt_weighted_max",
    "pocket_lddt_max",
    "pocket_lddt",
    "pocket_qcov_max",
    "pocket_qcov_weighted_max",
    "pocket_qcov",
]


@dataclass
class IngestConfig(_config.DataConfig):  # type: ignore
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
    """

    run_specific_stages: Any = ""
    skip_specific_stages: Any = ""

    def __post_init__(self) -> None:
        super().__post_init__()
        if isinstance(self.run_specific_stages, str):
            self.run_specific_stages = [
                stage for stage in self.run_specific_stages.split(",") if stage
            ]
        if isinstance(self.skip_specific_stages, str):
            self.skip_specific_stages = [
                stage for stage in self.skip_specific_stages.split(",") if stage
            ]


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
    score_ligand_level: bool = False
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
    artifact_within_entry_threshold: int = 0
    artifact_interacting_residue_count_threshold: int = 15
    # TODO-tjd: remove the two above


""" From
OleinikovasV
OleinikovasV commented Apr 22, 2024

Added updated artifacts list and curation, please, review the logic!

This no longer uses these at all:
artifact_within_entry_threshold: int = 15,
artifact_interacting_residue_count_threshold: int = 2,
"""


@dataclass
class ScatterConfig:
    """
    A class to represent batching parameters used to scatter data pipeline tasks.

    Attributes
    ----------
    two_char_batch_size: int
        Target number of two_char_codes per task batch.
    two_char_codes: str, default=None
        Only consider particular two character codes.
    pdb_ids : str, default=None
        Only consider particular pdb IDs.
    annotation_batch_size : int, default=220
        How many system annotations to generate in a given chunk
    skip_existing_entries : bool, default=True
        if the entry JSON already exists, skip generation
    make_entries_cpu : int, default=1
        misguided experiments in multiprocessing over C++ libs (bad idea)
    """

    two_char_batch_size: int = 2  # ~1060 codes / 500 workers
    two_char_codes: Any = ""
    pdb_ids: Any = ""
    annotation_batch_size: int = 220  # ~220000 PDB IDs / 1000 workers
    scorer_batch_size: int = 220
    skip_existing_entries: bool = True
    skip_missing_annotations: bool = True
    wipe_entries: bool = False
    wipe_annotations: bool = False
    make_entries_cpu: int = 4
    make_sub_dbs_cpu: int = 4
    make_scorers_cpu: int = 4
    download_alternative_datasets_threads: int = 10
    make_ligands_batch_size: int = 100
    number_id_col: str = "number_id_by_inchikeys"
    cluster_thresholds: list[int] = field(default_factory=lambda: [50, 70, 95, 100])
    skip_existing_clusters: bool = True
    stop_on_cluster: int = 0
    run_batch_searches_batch_size: int = 10050
    make_batch_scores_batch_size: int = 90
    split_config_dir: str = ""
    test_leakage: bool = False

    def __post_init__(self) -> None:
        for attr in [
            "two_char_codes",
            "pdb_ids",
        ]:
            value = getattr(self, attr, None)
            if isinstance(value, str):
                setattr(self, attr, [val for val in value.split(",") if val])


@dataclass
class LigandConfig:
    radius: int = 2
    nbits: int = 1024
    ligand_id_split_char: str = "__"
    save_top_k_similar_ligands: int = 5000
    multiply_by: int = 100
    score_name: str = "tanimoto_similarity_max"


SCHEMA = {
    "ingest": IngestConfig,
    "foldseek": FoldseekConfig,
    "mmseqs": MMSeqsConfig,
    "graph": GraphConfig,
    "annotation": AnnotationConfig,
    "entry": EntryConfig,
    "scorer": ScorerConfig,
    "scatter": ScatterConfig,
    "ligand": LigandConfig,
}


get_config = partial(_config._config, schema=SCHEMA, package_schema="data")
