# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0

from __future__ import annotations

from metaflow import FlowSpec, Parameter, kubernetes, environment, step, retry


MOUNT = "/plinder"
K8S = dict(
    cpu=1,
    image="us-east1-docker.pkg.dev/vantai-analysis/metaflow/plinder:v0.1.3-41-gd139d52b-dirty",
    node_selector={
        "topology.kubernetes.io/zone": "us-east1-b",
    },
    persistent_volume_claims={
        "plinder-data-gen-pvc": MOUNT,
    },
)
ENV = dict(
    vars=dict(
        PLINDER_MOUNT=MOUNT,
        PLINDER_RELEASE="2024-06",
        PLINDER_ITERATION="",
    )
)
LARGE_MEM = dict(
    cpu=7,
    memory=360000,
    tolerations=[
        dict(
            effect="NoSchedule",
            key="machine_type",
            value="n1-custom-380",
        )
    ]
)
WORKSTATION = dict(cpu=14, memory=14000)
PROTEIN_LEAKAGE = dict(
    cpu=10,
    memory=40000,
)
LIGAND_LEAKAGE = dict(
    cpu=4,
    memory=10000,
)


class PlinderSplitEvalFlow(FlowSpec):
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
        print(f"started split eval run with config: {self.config_file}")
        contents = gcs.download_as_str(gcs_path=self.config_file, bucket_name="plinder-collab-bucket")
        self.pipeline = IngestPipeline(conf=get_config(config_contents=contents))
        self.next(self.make_mmp_index)

    @kubernetes(**{**K8S, **WORKSTATION})
    @environment(**ENV)
    @retry
    @step
    def make_mmp_index(self):
        self.pipeline.make_mmp_index()
        self.next(self.scatter_make_splits)

    @kubernetes(**K8S)
    @environment(**ENV)
    @retry
    @step
    def scatter_make_splits(self):
        self.chunks = self.pipeline.scatter_make_splits()
        self.next(self.make_splits, foreach="chunks")

    @kubernetes(**{**K8S, **LARGE_MEM})
    @environment(**ENV)
    @retry
    @step
    def make_splits(self):
        self.pipeline.make_splits(self.input)
        self.next(self.join_make_splits)

    @kubernetes(**K8S)
    @environment(**ENV)
    @retry
    @step
    def join_make_splits(self, inputs):
        self.pipeline = inputs[0].pipeline
        self.merge_artifacts(inputs, exclude=["chunks"])
        self.next(self.end)

    # @kubernetes(**K8S)
    # @environment(**ENV)
    # @retry
    # @step
    # def scatter_compute_ligand_leakage(self):
    #     self.chunks = self.pipeline.scatter_compute_ligand_leakage()
    #     self.next(self.compute_ligand_leakage, foreach="chunks")
    #
    # @kubernetes(**{**K8S, **LIGAND_LEAKAGE})
    # @environment(**ENV)
    # @retry
    # @step
    # def compute_ligand_leakage(self):
    #     self.pipeline.compute_ligand_leakage(self.input)
    #     self.next(self.join_compute_ligand_leakage)
    #
    # @kubernetes(**K8S)
    # @environment(**ENV)
    # @retry
    # @step
    # def join_compute_ligand_leakage(self, inputs):
    #     self.pipeline = inputs[0].pipeline
    #     self.merge_artifacts(inputs, exclude=["chunks"])
    #     self.next(self.scatter_compute_protein_leakage)
    #
    # @kubernetes(**K8S)
    # @environment(**ENV)
    # @retry
    # @step
    # def scatter_compute_protein_leakage(self):
    #     self.chunks = self.pipeline.scatter_compute_protein_leakage()
    #     self.next(self.compute_protein_leakage, foreach="chunks")
    #
    # @kubernetes(**{**K8S, **PROTEIN_LEAKAGE})
    # @environment(**ENV)
    # @retry
    # @step
    # def compute_protein_leakage(self):
    #     self.pipeline.compute_protein_leakage(self.input)
    #     self.next(self.join_compute_protein_leakage)
    #
    # @kubernetes(**K8S)
    # @environment(**ENV)
    # @retry
    # @step
    # def join_compute_protein_leakage(self, inputs):
    #     self.pipeline = inputs[0].pipeline
    #     self.merge_artifacts(inputs, exclude=["chunks", "run_leakage"])
    #     self.next(self.end)

    @kubernetes(**K8S)
    @environment(**ENV)
    @retry
    @step
    def end(self):
        print("finished split eval run with:")
        print(f"pipeline config: {self.config_file}")
        print(f"split config dir: {self.pipeline.cfg.scatter.split_config_dir}")


if __name__ == "__main__":
    PlinderSplitEvalFlow()
