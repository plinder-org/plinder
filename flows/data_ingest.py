# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
"""
The filestore instance name is: plinder-data-gen.

"""
from __future__ import annotations

from metaflow import FlowSpec, Parameter, kubernetes, environment, step, retry


MOUNT = "/plinder"
K8S = dict(
    cpu=1,
    image="us-east1-docker.pkg.dev/vantai-analysis/vantai-rnd-images/plinder:v0.0.42-67-g6a138b8a",
    node_selector={
        "topology.kubernetes.io/zone": "us-east1-b",
    },
    persistent_volume_claims={
        "plinder-data-gen-pvc": MOUNT,
    },
)
ENV = dict(
    vars=dict(
        PLINDER_MOUNT="",
        PLINDER_RELEASE="2024-06",
        PLINDER_ITERATION="",
    )
)
DATABASES = dict(cpu=90, memory=82000)
WORKSTATION = dict(cpu=14, memory=14000)
WORKSTATION_MEM = dict(cpu=5, memory=48000)
LARGE_MEM = dict(
    cpu=1.75,
    memory=95000,
    tolerations=[
        dict(
            effect="NoSchedule",
            key="machine_type",
            value="n1-custom-380",
        )
    ]
)


class PlinderDataIngestFlow(FlowSpec):
    config_file = Parameter("config_file", required=True)

    @kubernetes(**K8S)
    @environment(**ENV)
    @retry
    @step
    def start(self):
        from plinder.core.utils import gcs
        from plinder.data.pipeline.config import get_config
        from plinder.data.pipeline.pipeline import IngestPipeline

        assert isinstance(self.config_file, str)
        if not self.config_file.startswith("gs:"):
            raise ValueError("--config_file must be a gs:// path")
        print(f"started data ingest run with config: {self.config_file}")
        contents = gcs.download_as_str(gcs_path=self.config_file, bucket_name="plinder-collab-bucket")
        self.pipeline = IngestPipeline(conf=get_config(config_contents=contents))
        self.next(self.scatter_make_components_and_communities)
    #     self.next(self.scatter_download_rcsb_files)
    #
    # @kubernetes(**K8S)
    # @environment(**ENV)
    # @retry
    # @step
    # def scatter_download_rcsb_files(self):
    #     self.chunks = self.pipeline.scatter_download_rcsb_files()
    #     self.next(self.download_rcsb_files, foreach="chunks")
    #
    # @kubernetes(**K8S)
    # @environment(**ENV)
    # @retry
    # @step
    # def download_rcsb_files(self):
    #     self.pipeline.download_rcsb_files(self.input)
    #     self.next(self.join_download_rcsb_files)
    #
    # @kubernetes(**K8S)
    # @environment(**ENV)
    # @retry
    # @step
    # def join_download_rcsb_files(self, inputs):
    #     self.merge_artifacts(inputs, exclude=["chunks"])
    #     self.next(self.download_alternative_datasets)
    #
    # @kubernetes(**K8S)
    # @environment(**ENV)
    # @retry
    # @step
    # def download_alternative_datasets(self):
    #     self.pipeline.download_alternative_datasets()
    #     self.next(self.make_dbs)
    #
    # @kubernetes(**{**K8S, **DATABASES})
    # @environment(**ENV)
    # @retry
    # @step
    # def make_dbs(self):
    #     self.pipeline.make_dbs()
    #     self.next(self.scatter_make_entries)
    #
    # @kubernetes(**K8S)
    # @environment(**ENV)
    # @retry
    # @step
    # def scatter_make_entries(self):
    #     self.chunks = self.pipeline.scatter_make_entries()
    #     self.next(self.make_entries, foreach="chunks")
    #
    # @kubernetes(**{**K8S, **WORKSTATION})
    # @environment(**ENV)
    # @retry
    # @step
    # def make_entries(self):
    #     self.pipeline.cfg.scatter.make_entries_cpu = WORKSTATION["cpu"]
    #     self.reruns = self.pipeline.make_entries(self.input)
    #     self.next(self.join_make_entries)
    #
    # @kubernetes(**K8S)
    # @environment(**ENV)
    # @retry
    # @step
    # def join_make_entries(self, inputs):
    #     self.merge_artifacts(inputs, exclude=["chunks", "reruns"])
    #     self.rerun = self.pipeline.join_make_entries([inp.reruns for inp in inputs])
    #     self.next(self.scatter_make_entries_second_try)
    #
    # @kubernetes(**K8S)
    # @environment(**ENV)
    # @retry
    # @step
    # def scatter_make_entries_second_try(self):
    #     self.original_pdb_ids = self.pipeline.cfg.scatter.pdb_ids
    #     self.pipeline.cfg.scatter.pdb_ids = self.rerun
    #     self.chunks = self.pipeline.scatter_make_entries()
    #     self.next(self.make_entries_second_try, foreach="chunks")
    #
    # @kubernetes(**{**K8S, **WORKSTATION_MEM})
    # @environment(**ENV)
    # @retry
    # @step
    # def make_entries_second_try(self):
    #     self.pipeline.cfg.scatter.make_entries_cpu = WORKSTATION_MEM["cpu"]
    #     self.pipeline.make_entries(self.input)
    #     self.next(self.join_make_entries_second_try)
    #
    # @kubernetes(**K8S)
    # @environment(**ENV)
    # @retry
    # @step
    # def join_make_entries_second_try(self, inputs):
    #     self.merge_artifacts(inputs, exclude=["chunks", "reruns"])
    #     self.pipeline.cfg.scatter.pdb_ids = self.original_pdb_ids
    #     self.next(self.scatter_structure_qc)
    #
    # @kubernetes(**K8S)
    # @environment(**ENV)
    # @retry
    # @step
    # def scatter_structure_qc(self):
    #     self.chunks = self.pipeline.scatter_structure_qc()
    #     self.next(self.structure_qc, foreach="chunks")
    #
    # @kubernetes(**{**K8S, **{"memory": 4000}})
    # @environment(**ENV)
    # @retry
    # @step
    # def structure_qc(self):
    #     self.pipeline.structure_qc(self.input)
    #     self.next(self.join_structure_qc)
    #
    # @kubernetes(**K8S)
    # @environment(**ENV)
    # @retry
    # @step
    # def join_structure_qc(self, inputs):
    #     self.merge_artifacts(inputs, exclude=["chunks"])
    #     self.next(self.scatter_make_system_archives)
    #
    # @kubernetes(**K8S)
    # @environment(**ENV)
    # @retry
    # @step
    # def scatter_make_system_archives(self):
    #     self.chunks = self.pipeline.scatter_make_system_archives()
    #     self.next(self.make_system_archives, foreach="chunks")
    #
    # @kubernetes(**K8S)
    # @environment(**ENV)
    # @retry
    # @step
    # def make_system_archives(self):
    #     self.pipeline.make_system_archives(self.input)
    #     self.next(self.join_make_system_archives)
    #
    # @kubernetes(**K8S)
    # @environment(**ENV)
    # @retry
    # @step
    # def join_make_system_archives(self, inputs):
    #     self.merge_artifacts(inputs, exclude=["chunks"])
    #     self.next(self.scatter_make_ligands)
    #
    # @kubernetes(**K8S)
    # @environment(**ENV)
    # @retry
    # @step
    # def scatter_make_ligands(self):
    #     self.chunks = self.pipeline.scatter_make_ligands()
    #     self.next(self.make_ligands, foreach="chunks")
    #
    # @kubernetes(**{**K8S, **WORKSTATION_MEM})
    # @environment(**ENV)
    # @retry
    # @step
    # def make_ligands(self):
    #     self.pipeline.make_ligands(self.input)
    #     self.next(self.join_make_ligands)
    #
    # @kubernetes(**K8S)
    # @environment(**ENV)
    # @retry
    # @step
    # def join_make_ligands(self, inputs):
    #     self.merge_artifacts(inputs, exclude=["chunks"])
    #     self.next(self.compute_ligand_fingerprints)
    #
    # @kubernetes(**{**K8S, **DATABASES})
    # @environment(**ENV)
    # @retry
    # @step
    # def compute_ligand_fingerprints(self):
    #     self.pipeline.compute_ligand_fingerprints()
    #     self.next(self.scatter_make_ligand_scores)
    #
    # @kubernetes(**K8S)
    # @environment(**ENV)
    # @retry
    # @step
    # def scatter_make_ligand_scores(self):
    #     self.chunks = self.pipeline.scatter_make_ligand_scores()
    #     self.next(self.make_ligand_scores, foreach="chunks")
    #
    # @kubernetes(**{**K8S, **WORKSTATION_MEM})
    # @environment(**ENV)
    # @retry
    # @step
    # def make_ligand_scores(self):
    #     self.pipeline.make_ligand_scores(self.input)
    #     self.reruns = self.pipeline.make_ligand_scores(self.input)
    #     self.next(self.join_make_ligand_scores)
    #
    # @kubernetes(**K8S)
    # @retry
    # @step
    # def join_make_ligand_scores(self, inputs):
    #     self.merge_artifacts(inputs, exclude=["chunks"])
    #     self.next(self.make_sub_dbs)
    #
    # @kubernetes(**{**K8S, **DATABASES})
    # @environment(**ENV)
    # @retry
    # @step
    # def make_sub_dbs(self):
    #     self.pipeline.make_sub_dbs()
    #     self.next(self.scatter_run_batch_searches)
    #
    # @kubernetes(**K8S)
    # @environment(**ENV)
    # @retry
    # @step
    # def scatter_run_batch_searches(self):
    #     self.chunks = self.pipeline.scatter_run_batch_searches()
    #     self.next(self.run_batch_searches, foreach="chunks")
    #
    # @kubernetes(**{**K8S, **DATABASES})
    # @environment(**ENV)
    # @retry
    # @step
    # def run_batch_searches(self):
    #     self.pipeline.run_batch_searches(self.input)
    #     self.next(self.join_run_batch_searches)
    #
    # @kubernetes(**K8S)
    # @environment(**ENV)
    # @retry
    # @step
    # def join_run_batch_searches(self, inputs):
    #     self.merge_artifacts(inputs, exclude=["chunks"])
    #     self.next(self.scatter_make_batch_scores)
    #
    # @kubernetes(**K8S)
    # @environment(**ENV)
    # @retry
    # @step
    # def scatter_make_batch_scores(self):
    #     self.chunks = self.pipeline.scatter_make_batch_scores()
    #     self.next(self.make_batch_scores, foreach="chunks")
    #
    # @kubernetes(**{**K8S, **{"memory": 10000}})
    # @environment(**ENV)
    # @retry
    # @step
    # def make_batch_scores(self):
    #     self.chunks = self.pipeline.make_batch_scores(self.input)
    #     self.next(self.join_make_batch_scores)
    #
    # @kubernetes(**K8S)
    # @environment(**ENV)
    # @retry
    # @step
    # def join_make_batch_scores(self, inputs):
    #     self.merge_artifacts(inputs, exclude=["chunks"])
    #     self.next(self.scatter_make_components_and_communities)

    @kubernetes(**K8S)
    @environment(**ENV)
    @retry
    @step
    def scatter_make_components_and_communities(self):
        self.chunks = self.pipeline.scatter_make_components_and_communities()
        self.next(self.make_components_and_communities, foreach="chunks")

    @kubernetes(**{**K8S, **LARGE_MEM})
    @environment(**ENV)
    @retry
    @step
    def make_components_and_communities(self):
        self.pipeline.make_components_and_communities(self.input)
        self.next(self.join_make_components_and_communities)

    @kubernetes(**K8S)
    @environment(**ENV)
    @retry
    @step
    def join_make_components_and_communities(self, inputs):
        self.pipeline = inputs[0].pipeline
        self.merge_artifacts(inputs, exclude=["chunks"])
        self.next(self.make_mmp_index)

    @kubernetes(**{**K8S, **WORKSTATION})
    @environment(**ENV)
    @retry
    @step
    def make_mmp_index(self):
        self.pipeline.make_mmp_index()
        self.next(self.end)

    @kubernetes(**K8S)
    @environment(**ENV)
    @retry
    @step
    def end(self):
        print(f"finished data ingest run with config: {self.config_file}")


if __name__ == "__main__":
    PlinderDataIngestFlow()
