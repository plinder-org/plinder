# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
"""
The filestore instance name is: plinder-data-gen.

TODO: The Metaflow pipeline still contains outdated V2 paths and has not been
tested end-to-end for V3.

"""
from __future__ import annotations

from metaflow import FlowSpec, Parameter, environment, kubernetes, retry, step

MOUNT = "/plinder"
K8S = dict(
    cpu=1,
    image="us-east1-docker.pkg.dev/vantai-analysis/metaflow/plinder:v0.2.2-63-g71bd2d22",
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
        PLINDER_RELEASE="2026-07",
        PLINDER_RELEASE_NUMBER="1",
    )
)
DATABASES = dict(cpu=90, memory=82000)
WORKSTATION = dict(cpu=14, memory=14000)
WORKSTATION_MEM = dict(cpu=5, memory=48000)
LARGE_MEM = dict(
    cpu=7,
    memory=380000,
    tolerations=[
        dict(
            effect="NoSchedule",
            key="machine_type",
            value="n1-custom-380",
        )
    ],
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
        contents = gcs.download_as_str(
            gcs_path=self.config_file, bucket_name="plinder-collab-bucket"
        )
        self.pipeline = IngestPipeline(conf=get_config(config_contents=contents))
        self.next(self.scatter_make_entries)

    @kubernetes(**{**K8S, **DATABASES})
    @environment(**ENV)
    @retry
    @step
    def make_dbs(self):
        self.pipeline.make_dbs()
        self.next(self.scatter_make_canonical_ligand_archives)

    @kubernetes(**K8S)
    @environment(**ENV)
    @retry
    @step
    def scatter_make_entries(self):
        self.chunks = self.pipeline.scatter_make_entries()
        self.next(self.make_entries, foreach="chunks")

    @kubernetes(**{**K8S, **WORKSTATION})
    @environment(**ENV)
    @retry
    @step
    def make_entries(self):
        self.pipeline.cfg.flow.make_entries_cpu = WORKSTATION["cpu"]
        self.reruns = self.pipeline.make_entries(self.input)
        self.next(self.join_make_entries)

    @kubernetes(**K8S)
    @environment(**ENV)
    @retry
    @step
    def join_make_entries(self, inputs):
        self.pipeline = inputs[0].pipeline
        self.merge_artifacts(inputs, exclude=["chunks", "reruns"])
        self.reruns = self.pipeline.join_make_entries(
            [input_.reruns for input_ in inputs]
        )
        self.next(self.scatter_collate_entries)

    @kubernetes(**K8S)
    @environment(**ENV)
    @retry
    @step
    def scatter_collate_entries(self):
        self.chunks = self.pipeline.scatter_collate_entries()
        self.next(self.collate_entries, foreach="chunks")

    @kubernetes(**{**K8S, **{"cpu": 2, "memory": 8000}})
    @environment(**ENV)
    @retry
    @step
    def collate_entries(self):
        self.pipeline.collate_entries(self.input)
        self.next(self.join_collate_entries)

    @kubernetes(**{**K8S, **WORKSTATION_MEM})
    @environment(**ENV)
    @retry
    @step
    def join_collate_entries(self, inputs):
        self.pipeline = inputs[0].pipeline
        self.merge_artifacts(inputs, exclude=["chunks"])
        self.pipeline.join_collate_entries([None for _ in inputs])
        self.next(self.make_protein_sequence_clusters)

    @kubernetes(**{**K8S, **DATABASES})
    @environment(**ENV)
    @retry
    @step
    def make_protein_sequence_clusters(self):
        self.pipeline.make_protein_sequence_clusters()
        self.next(self.make_dbs)

    @kubernetes(**K8S)
    @environment(**ENV)
    @retry
    @step
    def scatter_make_canonical_ligand_archives(self):
        self.chunks = self.pipeline.scatter_make_canonical_ligand_archives()
        self.next(self.make_canonical_ligand_archives, foreach="chunks")

    @kubernetes(**{**K8S, **{"memory": 4000}})
    @environment(**ENV)
    @retry
    @step
    def make_canonical_ligand_archives(self):
        self.pipeline.make_canonical_ligand_archives(self.input)
        self.next(self.join_make_canonical_ligand_archives)

    @kubernetes(**K8S)
    @environment(**ENV)
    @retry
    @step
    def join_make_canonical_ligand_archives(self, inputs):
        self.pipeline = inputs[0].pipeline
        self.merge_artifacts(inputs, exclude=["chunks"])
        self.next(self.finalize_ligand_archives)

    @kubernetes(**{**K8S, **WORKSTATION_MEM})
    @environment(**ENV)
    @retry
    @step
    def finalize_ligand_archives(self):
        self.pipeline.finalize_ligand_archives()
        self.next(self.compute_ligand_fingerprints)

    @kubernetes(**{**K8S, **DATABASES})
    @environment(**ENV)
    @retry
    @step
    def compute_ligand_fingerprints(self):
        self.pipeline.compute_ligand_fingerprints()
        self.next(self.scatter_make_ligand_scores)

    @kubernetes(**K8S)
    @environment(**ENV)
    @retry
    @step
    def scatter_make_ligand_scores(self):
        self.chunks = self.pipeline.scatter_make_ligand_scores()
        self.next(self.make_ligand_scores, foreach="chunks")

    @kubernetes(**{**K8S, **WORKSTATION_MEM})
    @environment(**ENV)
    @retry
    @step
    def make_ligand_scores(self):
        self.pipeline.make_ligand_scores(self.input)
        self.next(self.join_make_ligand_scores)

    @kubernetes(**K8S)
    @environment(**ENV)
    @retry
    @step
    def join_make_ligand_scores(self, inputs):
        self.pipeline = inputs[0].pipeline
        self.merge_artifacts(inputs, exclude=["chunks"])
        self.next(self.scatter_make_mhfp6_scores)

    @kubernetes(**K8S)
    @environment(**ENV)
    @retry
    @step
    def scatter_make_mhfp6_scores(self):
        self.chunks = self.pipeline.scatter_make_mhfp6_scores()
        self.next(self.make_mhfp6_scores, foreach="chunks")

    @kubernetes(**{**K8S, **WORKSTATION_MEM})
    @environment(**ENV)
    @retry
    @step
    def make_mhfp6_scores(self):
        self.pipeline.make_mhfp6_scores(self.input)
        self.next(self.join_make_mhfp6_scores)

    @kubernetes(**K8S)
    @environment(**ENV)
    @retry
    @step
    def join_make_mhfp6_scores(self, inputs):
        self.pipeline = inputs[0].pipeline
        self.merge_artifacts(inputs, exclude=["chunks"])
        self.next(self.annotate_ligand_similarity)

    @kubernetes(**{**K8S, **DATABASES})
    @environment(**ENV)
    @retry
    @step
    def annotate_ligand_similarity(self):
        self.pipeline.annotate_ligand_similarity()
        self.next(self.make_ligand_mmp_pairs)

    @kubernetes(**{**K8S, **WORKSTATION})
    @environment(**ENV)
    @retry
    @step
    def make_ligand_mmp_pairs(self):
        self.pipeline.make_ligand_mmp_pairs(threads=WORKSTATION["cpu"])
        self.next(self.make_sub_dbs)

    @kubernetes(**{**K8S, **DATABASES})
    @environment(**ENV)
    @retry
    @step
    def make_sub_dbs(self):
        self.pipeline.make_sub_dbs()
        self.next(self.scatter_run_batch_searches)

    @kubernetes(**K8S)
    @environment(**ENV)
    @retry
    @step
    def scatter_run_batch_searches(self):
        self.chunks = self.pipeline.scatter_run_batch_searches()
        self.next(self.run_batch_searches, foreach="chunks")

    @kubernetes(**{**K8S, **DATABASES})
    @environment(**ENV)
    @retry
    @step
    def run_batch_searches(self):
        self.pipeline.run_batch_searches(self.input)
        self.next(self.join_run_batch_searches)

    @kubernetes(**K8S)
    @environment(**ENV)
    @retry
    @step
    def join_run_batch_searches(self, inputs):
        self.pipeline = inputs[0].pipeline
        self.merge_artifacts(inputs, exclude=["chunks"])
        self.next(self.scatter_map_batch_alignments)

    @kubernetes(**K8S)
    @environment(**ENV)
    @retry
    @step
    def scatter_map_batch_alignments(self):
        self.chunks = self.pipeline.scatter_map_batch_alignments()
        self.next(self.map_batch_alignments, foreach="chunks")

    @kubernetes(**{**K8S, **{"memory": 10000}})
    @environment(**ENV)
    @retry
    @step
    def map_batch_alignments(self):
        self.pipeline.map_batch_alignments(self.input)
        self.next(self.join_map_batch_alignments)

    @kubernetes(**K8S)
    @environment(**ENV)
    @retry
    @step
    def join_map_batch_alignments(self, inputs):
        self.pipeline = inputs[0].pipeline
        self.merge_artifacts(inputs, exclude=["chunks"])
        self.next(self.scatter_collate_alignments)

    @kubernetes(**K8S)
    @environment(**ENV)
    @retry
    @step
    def scatter_collate_alignments(self):
        self.chunks = self.pipeline.scatter_collate_alignments()
        self.next(self.collate_alignments, foreach="chunks")

    @kubernetes(**{**K8S, **LARGE_MEM})
    @environment(**ENV)
    @retry
    @step
    def collate_alignments(self):
        self.pipeline.collate_alignments(self.input)
        self.next(self.join_collate_alignments)

    @kubernetes(**K8S)
    @environment(**ENV)
    @retry
    @step
    def join_collate_alignments(self, inputs):
        self.pipeline = inputs[0].pipeline
        self.merge_artifacts(inputs, exclude=["chunks"])
        self.next(self.finalize_alignments)

    @kubernetes(**{**K8S, **WORKSTATION_MEM})
    @environment(**ENV)
    @retry
    @step
    def finalize_alignments(self):
        self.pipeline.finalize_alignments()
        self.next(self.plan_score_batches)

    @kubernetes(**{**K8S, **WORKSTATION_MEM})
    @environment(**ENV)
    @retry
    @step
    def plan_score_batches(self):
        self.pipeline.plan_score_batches()
        self.next(self.plan_interface_scores)

    @kubernetes(**{**K8S, **WORKSTATION_MEM})
    @environment(**ENV)
    @retry
    @step
    def plan_interface_scores(self):
        self.pipeline.plan_interface_scores()
        self.next(self.scatter_make_interface_scores)

    @kubernetes(**K8S)
    @environment(**ENV)
    @retry
    @step
    def scatter_make_interface_scores(self):
        self.chunks = self.pipeline.scatter_make_interface_scores()
        self.next(self.make_interface_scores, foreach="chunks")

    @kubernetes(**{**K8S, **{"cpu": 4, "memory": 32000}})
    @environment(**ENV)
    @retry
    @step
    def make_interface_scores(self):
        self.pipeline.make_interface_scores(self.input)
        self.next(self.join_make_interface_scores)

    @kubernetes(**K8S)
    @environment(**ENV)
    @retry
    @step
    def join_make_interface_scores(self, inputs):
        self.pipeline = inputs[0].pipeline
        self.merge_artifacts(inputs, exclude=["chunks"])
        self.next(self.finalize_interface_scores)

    @kubernetes(**{**K8S, **WORKSTATION_MEM})
    @environment(**ENV)
    @retry
    @step
    def finalize_interface_scores(self):
        self.pipeline.finalize_interface_scores()
        self.next(self.scatter_make_batch_scores)

    @kubernetes(**K8S)
    @environment(**ENV)
    @retry
    @step
    def scatter_make_batch_scores(self):
        self.chunks = self.pipeline.scatter_make_batch_scores()
        self.next(self.make_batch_scores, foreach="chunks")

    @kubernetes(**{**K8S, **{"cpu": 4, "memory": 32000}})
    @environment(**ENV)
    @retry
    @step
    def make_batch_scores(self):
        self.pipeline.make_batch_scores(self.input)
        self.next(self.join_make_batch_scores)

    @kubernetes(**K8S)
    @environment(**ENV)
    @retry
    @step
    def join_make_batch_scores(self, inputs):
        self.pipeline = inputs[0].pipeline
        self.merge_artifacts(inputs, exclude=["chunks"])
        self.next(self.scatter_collate_ligand_3d_candidates)

    @kubernetes(**K8S)
    @environment(**ENV)
    @retry
    @step
    def scatter_collate_ligand_3d_candidates(self):
        self.chunks = self.pipeline.scatter_collate_ligand_3d_candidates()
        self.next(self.collate_ligand_3d_candidates, foreach="chunks")

    @kubernetes(**{**K8S, **{"cpu": 4, "memory": 32000}})
    @environment(**ENV)
    @retry
    @step
    def collate_ligand_3d_candidates(self):
        self.pipeline.collate_ligand_3d_candidates(self.input)
        self.next(self.join_collate_ligand_3d_candidates)

    @kubernetes(**K8S)
    @environment(**ENV)
    @retry
    @step
    def join_collate_ligand_3d_candidates(self, inputs):
        self.pipeline = inputs[0].pipeline
        self.merge_artifacts(inputs, exclude=["chunks"])
        self.next(self.plan_ligand_3d_scores)

    @kubernetes(**{**K8S, **WORKSTATION_MEM})
    @environment(**ENV)
    @retry
    @step
    def plan_ligand_3d_scores(self):
        self.pipeline.plan_ligand_3d_scores()
        self.next(self.scatter_make_ligand_3d_scores)

    @kubernetes(**K8S)
    @environment(**ENV)
    @retry
    @step
    def scatter_make_ligand_3d_scores(self):
        self.chunks = self.pipeline.scatter_make_ligand_3d_scores()
        self.next(self.make_ligand_3d_scores, foreach="chunks")

    @kubernetes(**{**K8S, **{"cpu": 4, "memory": 32000}})
    @environment(**ENV)
    @retry
    @step
    def make_ligand_3d_scores(self):
        self.pipeline.make_ligand_3d_scores(self.input)
        self.next(self.join_make_ligand_3d_scores)

    @kubernetes(**K8S)
    @environment(**ENV)
    @retry
    @step
    def join_make_ligand_3d_scores(self, inputs):
        self.pipeline = inputs[0].pipeline
        self.merge_artifacts(inputs, exclude=["chunks"])
        self.next(self.scatter_collate_ligand_3d_scores)

    @kubernetes(**K8S)
    @environment(**ENV)
    @retry
    @step
    def scatter_collate_ligand_3d_scores(self):
        self.chunks = self.pipeline.scatter_collate_ligand_3d_scores()
        self.next(self.collate_ligand_3d_scores, foreach="chunks")

    @kubernetes(**{**K8S, **{"cpu": 4, "memory": 32000}})
    @environment(**ENV)
    @retry
    @step
    def collate_ligand_3d_scores(self):
        self.pipeline.collate_ligand_3d_scores(self.input)
        self.next(self.join_collate_ligand_3d_scores)

    @kubernetes(**K8S)
    @environment(**ENV)
    @retry
    @step
    def join_collate_ligand_3d_scores(self, inputs):
        self.pipeline = inputs[0].pipeline
        self.merge_artifacts(inputs, exclude=["chunks"])
        self.next(self.scatter_merge_ligand_3d_scores)

    @kubernetes(**K8S)
    @environment(**ENV)
    @retry
    @step
    def scatter_merge_ligand_3d_scores(self):
        self.chunks = self.pipeline.scatter_merge_ligand_3d_scores()
        self.next(self.merge_ligand_3d_scores, foreach="chunks")

    @kubernetes(**{**K8S, **{"cpu": 4, "memory": 32000}})
    @environment(**ENV)
    @retry
    @step
    def merge_ligand_3d_scores(self):
        self.pipeline.merge_ligand_3d_scores(self.input)
        self.next(self.join_merge_ligand_3d_scores)

    @kubernetes(**K8S)
    @environment(**ENV)
    @retry
    @step
    def join_merge_ligand_3d_scores(self, inputs):
        self.pipeline = inputs[0].pipeline
        self.merge_artifacts(inputs, exclude=["chunks"])
        self.next(self.finalize_scores)

    @kubernetes(**{**K8S, **{"cpu": 4, "memory": 32000}})
    @environment(**ENV)
    @retry
    @step
    def finalize_scores(self):
        self.pipeline.finalize_scores()
        self.next(self.scatter_export_ligand_similarity_scores)

    @kubernetes(**K8S)
    @environment(**ENV)
    @retry
    @step
    def scatter_export_ligand_similarity_scores(self):
        self.chunks = self.pipeline.scatter_export_ligand_similarity_scores()
        self.next(self.export_ligand_similarity_scores, foreach="chunks")

    @kubernetes(**{**K8S, **{"cpu": 4, "memory": 32000}})
    @environment(**ENV)
    @retry
    @step
    def export_ligand_similarity_scores(self):
        self.pipeline.export_ligand_similarity_scores(self.input)
        self.next(self.join_export_ligand_similarity_scores)

    @kubernetes(**K8S)
    @environment(**ENV)
    @retry
    @step
    def join_export_ligand_similarity_scores(self, inputs):
        self.pipeline = inputs[0].pipeline
        self.merge_artifacts(inputs, exclude=["chunks"])
        self.next(self.finalize_ligand_similarity_scores)

    @kubernetes(**{**K8S, **LARGE_MEM})
    @environment(**ENV)
    @retry
    @step
    def finalize_ligand_similarity_scores(self):
        self.pipeline.finalize_ligand_similarity_scores()
        self.next(self.scatter_collate_partitions)

    @kubernetes(**K8S)
    @environment(**ENV)
    @retry
    @step
    def scatter_collate_partitions(self):
        self.chunks = self.pipeline.scatter_collate_partitions()
        self.next(self.collate_partitions, foreach="chunks")

    @kubernetes(**{**K8S, **LARGE_MEM})
    @environment(**ENV)
    @retry
    @step
    def collate_partitions(self):
        self.pipeline.collate_partitions(self.input)
        self.next(self.join_collate_partitions)

    @kubernetes(**K8S)
    @environment(**ENV)
    @retry
    @step
    def join_collate_partitions(self, inputs):
        self.pipeline = inputs[0].pipeline
        self.merge_artifacts(inputs, exclude=["chunks"])
        self.next(self.make_linked_apo_structures)

    @kubernetes(**{**K8S, **LARGE_MEM})
    @environment(**ENV)
    @retry
    @step
    def make_linked_apo_structures(self):
        self.pipeline.make_linked_apo_structures()
        self.next(self.plan_clusters)

    @kubernetes(**K8S)
    @environment(**ENV)
    @retry
    @step
    def plan_clusters(self):
        self.pipeline.plan_clusters()
        self.next(self.scatter_make_symmetric_edge_fragments)

    @kubernetes(**K8S)
    @environment(**ENV)
    @retry
    @step
    def scatter_make_symmetric_edge_fragments(self):
        self.chunks = self.pipeline.scatter_make_symmetric_edge_fragments()
        self.next(self.make_symmetric_edge_fragments, foreach="chunks")

    @kubernetes(**{**K8S, **LARGE_MEM})
    @environment(**ENV)
    @retry
    @step
    def make_symmetric_edge_fragments(self):
        self.pipeline.make_symmetric_edge_fragments(self.input)
        self.next(self.join_make_symmetric_edge_fragments)

    @kubernetes(**K8S)
    @environment(**ENV)
    @retry
    @step
    def join_make_symmetric_edge_fragments(self, inputs):
        self.pipeline = inputs[0].pipeline
        self.merge_artifacts(inputs, exclude=["chunks"])
        self.next(self.scatter_make_symmetric_edge_shards)

    @kubernetes(**K8S)
    @environment(**ENV)
    @retry
    @step
    def scatter_make_symmetric_edge_shards(self):
        self.chunks = self.pipeline.scatter_make_symmetric_edge_shards()
        self.next(self.make_symmetric_edge_shards, foreach="chunks")

    @kubernetes(**{**K8S, **LARGE_MEM})
    @environment(**ENV)
    @retry
    @step
    def make_symmetric_edge_shards(self):
        self.pipeline.make_symmetric_edge_shards(self.input)
        self.next(self.join_make_symmetric_edge_shards)

    @kubernetes(**K8S)
    @environment(**ENV)
    @retry
    @step
    def join_make_symmetric_edge_shards(self, inputs):
        self.pipeline = inputs[0].pipeline
        self.merge_artifacts(inputs, exclude=["chunks"])
        self.next(self.scatter_make_component_reductions)

    @kubernetes(**K8S)
    @environment(**ENV)
    @retry
    @step
    def scatter_make_component_reductions(self):
        self.chunks = self.pipeline.scatter_make_component_reductions()
        self.next(self.make_component_reductions, foreach="chunks")

    @kubernetes(**{**K8S, **LARGE_MEM})
    @environment(**ENV)
    @retry
    @step
    def make_component_reductions(self):
        self.pipeline.make_component_reductions(self.input)
        self.next(self.join_make_component_reductions)

    @kubernetes(**{**K8S, **LARGE_MEM})
    @environment(**ENV)
    @retry
    @step
    def join_make_component_reductions(self, inputs):
        self.pipeline = inputs[0].pipeline
        self.merge_artifacts(inputs, exclude=["chunks"])
        self.pipeline.merge_component_reductions()
        self.next(self.scatter_make_set_covers)

    @kubernetes(**K8S)
    @environment(**ENV)
    @retry
    @step
    def scatter_make_set_covers(self):
        self.chunks = self.pipeline.scatter_make_set_covers()
        self.next(self.make_set_covers, foreach="chunks")

    @kubernetes(**{**K8S, **LARGE_MEM})
    @environment(**ENV)
    @retry
    @step
    def make_set_covers(self):
        self.pipeline.make_set_covers(self.input)
        self.next(self.join_make_set_covers)

    @kubernetes(**K8S)
    @environment(**ENV)
    @retry
    @step
    def join_make_set_covers(self, inputs):
        self.pipeline = inputs[0].pipeline
        self.merge_artifacts(inputs, exclude=["chunks"])
        self.next(self.scatter_make_directed_set_covers)

    @kubernetes(**K8S)
    @environment(**ENV)
    @retry
    @step
    def scatter_make_directed_set_covers(self):
        self.chunks = self.pipeline.scatter_make_directed_set_covers()
        self.next(self.make_directed_set_covers, foreach="chunks")

    @kubernetes(**{**K8S, **LARGE_MEM})
    @environment(**ENV)
    @retry
    @step
    def make_directed_set_covers(self):
        self.pipeline.make_directed_set_covers(self.input)
        self.next(self.join_make_directed_set_covers)

    @kubernetes(**K8S)
    @environment(**ENV)
    @retry
    @step
    def join_make_directed_set_covers(self, inputs):
        self.pipeline = inputs[0].pipeline
        self.merge_artifacts(inputs, exclude=["chunks"])
        self.next(self.summarize_clusters)

    @kubernetes(**{**K8S, **LARGE_MEM})
    @environment(**ENV)
    @retry
    @step
    def summarize_clusters(self):
        self.pipeline.summarize_clusters()
        self.next(self.finalize_index)

    @kubernetes(**{**K8S, **LARGE_MEM})
    @environment(**ENV)
    @retry
    @step
    def finalize_index(self):
        self.pipeline.finalize_index()
        self.next(self.end)

    @kubernetes(**K8S)
    @environment(**ENV)
    @retry
    @step
    def end(self):
        print(f"finished data ingest run with config: {self.config_file}")


if __name__ == "__main__":
    PlinderDataIngestFlow()
