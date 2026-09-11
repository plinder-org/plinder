# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
import tempfile
from pathlib import Path
from typing import Any, Literal, Optional

import pandas as pd
from omegaconf import DictConfig

from plinder.core.utils.log import setup_logger
from plinder.data.pipeline import config, tasks, utils
from plinder.data.pipeline.ingest import resolve_source_roots

LOG = setup_logger(__name__)


class IngestPipeline:
    """
    Mimic the required metaflow DAG pattern of
        - scatter
        - compute
        - join
    by convention in method names.

    scatter methods return lists of lists of primitives
    compute methods may return something if intended to be joined
    join methods may return something if used elsewhere

    Note
    ----
    business logic is implemented in tasks.py. The
    Pipeline is merely an interface between configuration
    and functions, which can be mirrored in metaflow.
    """

    def __init__(
        self,
        conf: Optional[DictConfig] = None,
        config_file: Optional[str] = None,
        config_contents: Optional[str] = None,
        config_args: Optional[list[str]] = None,
        cached: bool = True,
    ):
        self.cfg = config.get_config(
            config=conf,
            config_file=config_file,
            config_contents=config_contents,
            config_args=config_args,
            cached=cached,
        )
        self.plinder_dir = Path(self.cfg.data.plinder_dir)
        LOG.info(f"plinder_dir={self.plinder_dir}")

    def _entry_source_roots(self) -> tuple[Path, Path]:
        """Resolve configured source archives for entry generation."""
        return resolve_source_roots(
            data_dir=self.plinder_dir,
            cif_root=self.cfg.source.pdb_nextgen_root or None,
            validation_root=self.cfg.source.validation_root or None,
        )

    def _seqres_source(self) -> Path | None:
        """Resolve an optional pre-fetched PDB SEQRES input."""
        configured = str(self.cfg.source.seqres_path)
        if not configured:
            return None
        path = Path(configured).expanduser()
        if not path.is_absolute():
            path = self.plinder_dir / path
        return path.resolve()

    def __setstate__(self, state: dict[str, Any]) -> None:
        self.cfg = config.get_config(config=state.pop("cfg"))
        for k, v in state.items():
            setattr(self, k, v)
        self.plinder_dir = Path(self.cfg.data.plinder_dir)
        LOG.info(f"plinder_dir={self.plinder_dir}")

    @utils.ingest_flow_control
    def scatter_download_rcsb_files(self) -> list[list[str]]:
        chunks: list[list[str]] = tasks.scatter_download_rcsb_files(
            data_dir=self.plinder_dir,
            two_char_codes=self.cfg.context.two_char_codes,
            batch_size=self.cfg.flow.download_rcsb_files_batch_size,
        )
        return chunks

    @utils.ingest_flow_control
    def download_rcsb_files(self, two_char_codes: list[str]) -> None:
        tasks.download_rcsb_files(
            data_dir=self.plinder_dir,
            two_char_codes=two_char_codes,
        )

    @utils.ingest_flow_control
    def download_alternative_datasets(self) -> None:
        tasks.download_alternative_datasets(
            data_dir=self.plinder_dir,
            threads=self.cfg.flow.download_alternative_datasets_threads,
            force_update=self.cfg.data.force_update,
        )

    @utils.ingest_flow_control
    def make_protein_sequence_clusters(self) -> None:
        tasks.make_protein_sequence_clusters(
            data_dir=self.plinder_dir,
            scratch_dir=Path(tempfile.gettempdir()),
            cpu=self.cfg.flow.protein_clustering_cpu,
            identity=self.cfg.flow.protein_sequence_cluster_identity,
            coverage=self.cfg.flow.protein_cluster_coverage,
            force_update=self.cfg.data.force_update,
        )

    @utils.ingest_flow_control
    def make_dbs(self) -> None:
        cif_root, _ = self._entry_source_roots()
        from plinder.data.pipeline.score import (
            make_foldseek_input_manifest,
            plan_protein_scoring,
        )

        if self.cfg.foldseek.max_seqs != self.cfg.mmseqs.max_seqs:
            raise ValueError(
                "Foldseek and MMseqs max_seqs must match in one scoring plan"
            )
        plan_protein_scoring(
            self.plinder_dir,
            pdb_ids=list(self.cfg.context.pdb_ids),
            two_char_codes=list(self.cfg.context.two_char_codes),
            max_seqs=self.cfg.foldseek.max_seqs,
        )
        cif_root = make_foldseek_input_manifest(self.plinder_dir, cif_root)
        scratch_root = Path(tempfile.gettempdir()) / "plinder-full-search-dbs"
        tasks.make_dbs(
            data_dir=self.plinder_dir,
            sub_databases=self.cfg.scorer.sub_databases,
            cpu=self.cfg.flow.make_dbs_cpu,
            cif_root=cif_root,
            seqres_path=self._seqres_source(),
            scratch_dir=(scratch_root / "work" if scratch_root is not None else None),
            build_dir=(scratch_root / "build" if scratch_root is not None else None),
            index=False,
            force_update=self.cfg.data.force_update,
        )

    @utils.ingest_flow_control
    def scatter_make_entries(self) -> list[list[str]]:
        force_update = (
            self.cfg.data.force_update or self.cfg.flow.make_entries_force_update
        )
        cif_root, validation_root = self._entry_source_roots()
        chunks: list[list[str]] = tasks.scatter_make_entries(
            data_dir=self.plinder_dir,
            cif_root=cif_root,
            validation_root=validation_root,
            batch_size=self.cfg.flow.make_entries_batch_size,
            two_char_codes=self.cfg.context.two_char_codes,
            pdb_ids=self.cfg.context.pdb_ids,
            force_update=force_update,
            discovery_threads=self.cfg.source.discovery_threads,
            interface_min_residues=self.cfg.interface.min_interface_residues,
            interface_annotate_prodigy=self.cfg.interface.annotate_prodigy,
            ingest_mode=self.cfg.flow.make_entries_mode,
        )
        return chunks

    @utils.ingest_flow_control
    def make_entries(self, pdb_ids: list[str]) -> list[str]:
        force_update = (
            self.cfg.data.force_update or self.cfg.flow.make_entries_force_update
        )
        cif_root, validation_root = self._entry_source_roots()
        failed: list[str] = tasks.make_entries(
            data_dir=self.plinder_dir,
            pdb_ids=pdb_ids,
            cif_root=cif_root,
            validation_root=validation_root,
            force_update=force_update,
            cpu=self.cfg.flow.make_entries_cpu,
            annotation_cfg=self.cfg.annotation,
            entry_cfg=self.cfg.entry,
            interface_cfg=self.cfg.interface,
            ingest_mode=self.cfg.flow.make_entries_mode,
        )
        return failed

    @utils.ingest_flow_control
    def join_make_entries(self, reruns: list[list[str]]) -> list[str]:
        catted = []
        for rerun in reruns:
            catted.extend([item[-4:] for item in rerun])
        return catted

    @utils.ingest_flow_control
    def scatter_collate_entries(self) -> list[list[str]]:
        return tasks.scatter_collate_entries(
            data_dir=self.plinder_dir,
            batch_size=self.cfg.flow.collate_entries_batch_size,
        )

    @utils.ingest_flow_control
    def collate_entries(self, two_char_codes: list[str]) -> None:
        tasks.collate_entries(
            data_dir=self.plinder_dir,
            two_char_codes=two_char_codes,
            cpu=self.cfg.flow.collate_entries_cpu,
            memory_limit=self.cfg.flow.collate_entries_memory_limit,
        )

    @utils.ingest_flow_control
    def join_collate_entries(self, outputs: list[None]) -> None:
        del outputs
        tasks.finalize_entry_collation(
            data_dir=self.plinder_dir,
            cpu=self.cfg.flow.finalize_entries_cpu,
            memory_limit=self.cfg.flow.finalize_entries_memory_limit,
        )

    @utils.ingest_flow_control
    def compute_ligand_fingerprints(self) -> None:
        tasks.compute_ligand_fingerprints(
            data_dir=self.plinder_dir,
            cofactor_similarity_threshold=self.cfg.ligand.cofactor_similarity_threshold,
        )

    @utils.ingest_flow_control
    def scatter_make_ligand_scores(self) -> list[list[int]]:
        ligand_ids: list[list[int]] = tasks.scatter_make_ligand_scores(
            data_dir=self.plinder_dir,
            batch_size=self.cfg.flow.make_ligands_batch_size,
            number_id_col=self.cfg.ligand.number_id_col,
        )
        return ligand_ids

    @utils.ingest_flow_control
    def make_ligand_scores(self, ligand_ids: list[int]) -> None:
        tasks.make_ligand_scores(
            data_dir=self.plinder_dir,
            ligand_ids=ligand_ids,
            minimum_similarity=self.cfg.ligand.minimum_similarity,
            number_id_col=self.cfg.ligand.number_id_col,
        )

    @utils.ingest_flow_control
    def scatter_make_mhfp6_scores(self) -> list[list[int]]:
        # MHFP6 clustering is opt-in: only score when the metric is configured,
        # so dropping it from cluster_metrics also skips its all-pairs scoring.
        from plinder.data.annotations.get_similarity_scores import MHFP6_METRIC

        if MHFP6_METRIC not in self.cfg.flow.cluster_metrics:
            LOG.info(
                "scatter_make_mhfp6_scores: MHFP6 not in cluster_metrics; skipping"
            )
            return [[]]
        ligand_ids: list[list[int]] = tasks.scatter_make_mhfp6_scores(
            data_dir=self.plinder_dir,
            batch_size=self.cfg.flow.make_ligands_batch_size,
            number_id_col=self.cfg.ligand.number_id_col,
        )
        return ligand_ids

    @utils.ingest_flow_control
    def make_mhfp6_scores(self, ligand_ids: list[int]) -> None:
        tasks.make_mhfp6_scores(
            data_dir=self.plinder_dir,
            ligand_ids=ligand_ids,
            minimum_similarity=self.cfg.ligand.minimum_similarity,
            number_id_col=self.cfg.ligand.number_id_col,
        )

    @utils.ingest_flow_control
    def annotate_ligand_similarity(self) -> None:
        tasks.annotate_ligand_similarity(data_dir=self.plinder_dir)

    @utils.ingest_flow_control
    def make_ligand_mmp_pairs(self, *, threads: int | None = None) -> None:
        tasks.make_ligand_mmp_pairs(
            data_dir=self.plinder_dir,
            scratch_dir=Path(tempfile.gettempdir()) / "plinder-mmp",
            threads=threads or self.cfg.flow.clustering_cpu,
            force_update=self.cfg.data.force_update,
        )

    @utils.ingest_flow_control
    def scatter_make_canonical_ligand_archives(self) -> list[list[str]]:
        chunks: list[list[str]] = tasks.scatter_make_canonical_ligand_archives(
            data_dir=self.plinder_dir,
            two_char_codes=self.cfg.context.two_char_codes,
            pdb_ids=self.cfg.context.pdb_ids,
            batch_size=self.cfg.flow.download_rcsb_files_batch_size,
        )
        return chunks

    @utils.ingest_flow_control
    def make_canonical_ligand_archives(self, two_char_codes: list[str]) -> None:
        tasks.make_canonical_ligand_archives(
            data_dir=self.plinder_dir,
            two_char_codes=two_char_codes,
        )

    @utils.ingest_flow_control
    def finalize_ligand_archives(self) -> None:
        from plinder.data.pipeline.score import finalize_ligand_archives

        finalize_ligand_archives(self.plinder_dir)

    @utils.ingest_flow_control
    def make_sub_dbs(self) -> None:
        tasks.make_sub_dbs(
            data_dir=self.plinder_dir,
            sub_databases=self.cfg.scorer.sub_databases,
            cpu=self.cfg.flow.make_sub_dbs_cpu,
            scratch_dir=Path(tempfile.gettempdir()) / "plinder-exact-search-dbs",
        )

    @utils.ingest_flow_control
    def scatter_run_batch_searches(self) -> list[list[str]]:
        chunks: list[list[str]] = tasks.scatter_protein_scoring(
            data_dir=self.plinder_dir,
            batch_size=self.cfg.flow.run_batch_searches_batch_size,
            two_char_codes=self.cfg.context.two_char_codes,
            pdb_ids=self.cfg.context.pdb_ids,
            search_dbs=self.cfg.scorer.sub_databases,
        )
        return chunks

    @utils.ingest_flow_control
    def run_batch_searches(self, pdb_ids: list[str]) -> None:
        tasks.run_batch_searches(
            data_dir=self.plinder_dir,
            pdb_ids=pdb_ids,
            scorer_cfg=self.cfg.scorer,
            foldseek_cfg=self.cfg.foldseek,
            mmseqs_cfg=self.cfg.mmseqs,
            cpu=self.cfg.flow.make_scorers_cpu,
            force_update=self.cfg.data.force_update,
        )

    @utils.ingest_flow_control
    def scatter_map_batch_alignments(self) -> list[tuple[str, list[str]]]:
        work: list[tuple[str, list[str]]] = []
        for search_db in self.cfg.scorer.sub_databases:
            chunks = tasks.scatter_missing_alignment_mappings(
                data_dir=self.plinder_dir,
                batch_size=self.cfg.flow.map_batch_alignments_batch_size,
                search_db=search_db,
            )
            work.extend((search_db, shards) for shards in chunks if shards)
        return work or [(str(self.cfg.scorer.sub_databases[0]), [])]

    @utils.ingest_flow_control
    def map_batch_alignments(self, work: tuple[str, list[str]]) -> None:
        search_db, shards = work
        force_update = self.cfg.data.force_update
        tasks.map_batch_alignments(
            data_dir=self.plinder_dir,
            shards=shards,
            scorer_cfg=self.cfg.scorer,
            force_update=force_update,
            search_db=search_db,
        )

    @utils.ingest_flow_control
    def scatter_make_batch_scores(self) -> list[list[str]]:
        chunks: list[list[str]] = tasks.scatter_missing_scores(
            data_dir=self.plinder_dir,
            batch_size=self.cfg.flow.make_batch_scores_batch_size,
            scorer_cfg=self.cfg.scorer,
            search_dbs=self.cfg.scorer.sub_databases,
        )
        return chunks

    @utils.ingest_flow_control
    def make_batch_scores(self, pdb_ids: list[str]) -> None:
        force_update = (
            self.cfg.data.force_update or self.cfg.flow.make_batch_scores_force_update
        )
        tasks.make_batch_scores(
            data_dir=self.plinder_dir,
            pdb_ids=pdb_ids,
            scorer_cfg=self.cfg.scorer,
            force_update=force_update,
            threads=self.cfg.flow.make_batch_scores_cpu,
            defer_ligand_3d=True,
        )

    @utils.ingest_flow_control
    def scatter_collate_ligand_3d_candidates(self) -> list[list[str]]:
        return tasks.scatter_ligand_3d_candidate_shards(
            data_dir=self.plinder_dir,
            batch_size=self.cfg.flow.collate_ligand_3d_candidates_batch_size,
        )

    @utils.ingest_flow_control
    def collate_ligand_3d_candidates(self, shards: list[str]) -> None:
        tasks.collate_ligand_3d_candidates(
            data_dir=self.plinder_dir,
            shards=shards,
            scratch_dir=Path(tempfile.gettempdir()) / "plinder-ligand-3d-candidates",
            threads=self.cfg.flow.make_ligand_3d_scores_cpu,
        )

    @utils.ingest_flow_control
    def plan_ligand_3d_scores(self) -> None:
        from plinder.data.pipeline.score import plan_ligand_3d_batches

        plan_ligand_3d_batches(
            self.plinder_dir,
            batch_size=self.cfg.flow.make_ligand_3d_scores_batch_size,
            threads=self.cfg.flow.make_ligand_3d_scores_cpu,
            scratch_dir=Path(tempfile.gettempdir()) / "plinder-ligand-3d-plan",
        )

    @utils.ingest_flow_control
    def scatter_make_ligand_3d_scores(self) -> list[list[int]]:
        return tasks.scatter_ligand_3d_score_batches(data_dir=self.plinder_dir)

    @utils.ingest_flow_control
    def make_ligand_3d_scores(self, batch_indices: list[int]) -> None:
        from plinder.data.pipeline.score import _ligand_3d_batch

        force_update = (
            self.cfg.data.force_update or self.cfg.flow.make_batch_scores_force_update
        )
        for batch_index in batch_indices:
            pairs = _ligand_3d_batch(
                self.plinder_dir,
                batch_index,
                self.cfg.flow.make_ligand_3d_scores_batch_size,
            )
            tasks.make_ligand_3d_scores(
                data_dir=self.plinder_dir,
                pairs=pairs,
                batch_index=batch_index,
                scorer_cfg=self.cfg.scorer,
                force_update=force_update,
                scratch_dir=(
                    Path(tempfile.gettempdir())
                    / "plinder-ligand-3d-score"
                    / str(batch_index)
                ),
                threads=self.cfg.flow.make_ligand_3d_scores_cpu,
            )

    @utils.ingest_flow_control
    def scatter_collate_ligand_3d_scores(self) -> list[list[str]]:
        return tasks.scatter_ligand_3d_query_shards(
            data_dir=self.plinder_dir,
            batch_size=self.cfg.flow.collate_ligand_3d_scores_batch_size,
        )

    @utils.ingest_flow_control
    def collate_ligand_3d_scores(self, shards: list[str]) -> None:
        tasks.collate_ligand_3d_scores(
            data_dir=self.plinder_dir,
            shards=shards,
            scratch_dir=Path(tempfile.gettempdir()) / "plinder-ligand-3d-collate",
            threads=self.cfg.flow.make_ligand_3d_scores_cpu,
        )

    @utils.ingest_flow_control
    def scatter_merge_ligand_3d_scores(self) -> list[list[str]]:
        return tasks.scatter_ligand_3d_merge(
            data_dir=self.plinder_dir,
            batch_size=self.cfg.flow.merge_ligand_3d_scores_batch_size,
        )

    @utils.ingest_flow_control
    def merge_ligand_3d_scores(self, shards: list[str]) -> None:
        force_update = (
            self.cfg.data.force_update or self.cfg.flow.make_batch_scores_force_update
        )
        tasks.merge_ligand_3d_scores(
            data_dir=self.plinder_dir,
            shards=shards,
            scorer_cfg=self.cfg.scorer,
            force_update=force_update,
            scratch_dir=Path(tempfile.gettempdir()) / "plinder-ligand-3d-merge",
            threads=self.cfg.flow.make_ligand_3d_scores_cpu,
        )

    @utils.ingest_flow_control
    def finalize_scores(self) -> None:
        from plinder.data.pipeline.score import finalize_ligand_3d_artifacts

        finalize_ligand_3d_artifacts(self.plinder_dir)

    @utils.ingest_flow_control
    def scatter_export_ligand_similarity_scores(self) -> list[list[str]]:
        return tasks.scatter_ligand_3d_query_shards(
            data_dir=self.plinder_dir,
            batch_size=self.cfg.flow.collate_ligand_3d_scores_batch_size,
        )

    @utils.ingest_flow_control
    def export_ligand_similarity_scores(self, shards: list[str]) -> None:
        from plinder.data.pipeline.score import export_ligand_similarity_scores_batch

        export_ligand_similarity_scores_batch(
            self.plinder_dir,
            output_dir=(self.plinder_dir / "exports" / "ligand_similarity_scores"),
            shards=shards,
            scratch_dir=(
                Path(tempfile.gettempdir()) / "plinder-ligand-similarity-export"
            ),
            threads=self.cfg.flow.make_ligand_3d_scores_cpu,
        )

    @utils.ingest_flow_control
    def finalize_ligand_similarity_scores(self) -> None:
        from plinder.data.pipeline.score import finalize_ligand_similarity_scores

        finalize_ligand_similarity_scores(
            self.plinder_dir,
            source_dir=(self.plinder_dir / "exports" / "ligand_similarity_scores"),
            output=self.plinder_dir / "exports" / "ligand_similarity_scores.parquet",
            scratch_dir=(
                Path(tempfile.gettempdir()) / "plinder-ligand-similarity-finalize"
            ),
            threads=self.cfg.flow.make_ligand_3d_scores_cpu,
        )

    @utils.ingest_flow_control
    def scatter_collate_alignments(self) -> list[list[str]]:
        return tasks.scatter_collate_alignments(data_dir=self.plinder_dir)

    @utils.ingest_flow_control
    def collate_alignments(self, partition: list[str]) -> None:
        tasks.collate_alignments(data_dir=self.plinder_dir, partition=partition)

    @utils.ingest_flow_control
    def finalize_alignments(self) -> None:
        from plinder.data.pipeline.score import finalize_alignment_artifacts

        finalize_alignment_artifacts(self.plinder_dir)

    @utils.ingest_flow_control
    def plan_score_batches(self) -> None:
        from plinder.data.pipeline.score import plan_score_batches

        plan_score_batches(
            self.plinder_dir,
            batch_size=self.cfg.flow.make_batch_scores_batch_size,
            threads=self.cfg.flow.make_batch_scores_cpu,
            scratch_dir=Path(tempfile.gettempdir()) / "plinder-score-plan",
            max_query_protein_chains=self.cfg.scorer.max_query_protein_chains,
            max_query_proper_ligand_chains=(
                self.cfg.scorer.max_query_proper_ligand_chains
            ),
        )

    @utils.ingest_flow_control
    def plan_interface_scores(self) -> None:
        from plinder.data.pipeline.score import plan_interface_scoring

        plan_interface_scoring(
            self.plinder_dir,
            batch_size=self.cfg.flow.make_interface_scores_batch_size,
        )

    @utils.ingest_flow_control
    def scatter_make_interface_scores(self) -> list[list[str]]:
        from plinder.data.pipeline.score import (
            INTERFACE_SCORE_WORK_RELATIVE,
            _load_interface_score_plan,
        )

        plan = _load_interface_score_plan(self.plinder_dir)
        work = (
            pd.read_parquet(
                self.plinder_dir / INTERFACE_SCORE_WORK_RELATIVE,
                columns=["shard"],
            )["shard"]
            .astype(str)
            .tolist()
        )
        batch_size = int(plan["batch_size"])
        return [
            work[start : start + batch_size]
            for start in range(0, len(work), batch_size)
        ] or [[]]

    @utils.ingest_flow_control
    def make_interface_scores(self, shards: list[str]) -> None:
        if not shards:
            return
        from plinder.data.pipeline.score import score_interface_qcov_shards

        score_interface_qcov_shards(
            self.plinder_dir,
            shards=shards,
            scratch_dir=Path(tempfile.gettempdir()) / "plinder-interface-scores",
            threads=self.cfg.flow.make_interface_scores_cpu,
            memory_limit=self.cfg.flow.make_interface_scores_memory_limit,
            force_update=self.cfg.data.force_update,
        )

    @utils.ingest_flow_control
    def finalize_interface_scores(self) -> None:
        from plinder.data.pipeline.score import finalize_interface_similarity_scores

        finalize_interface_similarity_scores(
            self.plinder_dir,
            scratch_dir=Path(tempfile.gettempdir()) / "plinder-interface-finalize",
            threads=self.cfg.flow.make_interface_scores_cpu,
            memory_limit=self.cfg.flow.make_interface_scores_memory_limit,
        )

    @utils.ingest_flow_control
    def scatter_collate_partitions(self) -> list[list[str]]:
        chunks: list[list[str]] = tasks.scatter_collate_partitions()
        # Holo scores are already published as two-character query shards by
        # merge_ligand_3d_scores. Only linked apo/pred scores still use this
        # partition collation path.
        chunks = [
            chunk
            for chunk in chunks
            if chunk
            and chunk[0] in self.cfg.scorer.sub_databases
            and chunk[0] in {"apo", "pred"}
        ]
        return chunks or [[]]

    @utils.ingest_flow_control
    def collate_partitions(self, partition: list[str]) -> None:
        tasks.collate_partitions(
            data_dir=self.plinder_dir,
            partition=partition,
            scratch_dir=Path(tempfile.gettempdir()) / "plinder-score-partitions",
            threads=self.cfg.flow.collate_partitions_cpu,
            memory_limit=self.cfg.flow.collate_partitions_memory_limit,
        )

    @utils.ingest_flow_control
    def make_linked_apo_structures(self) -> None:
        if "apo" not in self.cfg.scorer.sub_databases:
            LOG.info("make_linked_apo_structures: apo scoring is disabled")
            return
        tasks.make_linked_apo_structures(
            data_dir=self.plinder_dir,
            scratch_dir=Path(tempfile.gettempdir()) / "plinder-linked-apo",
            threads=self.cfg.flow.collate_partitions_cpu,
            memory_limit=self.cfg.flow.collate_partitions_memory_limit,
        )

    def _cluster_entities(
        self,
    ) -> list[tuple[Literal["ligand", "interface"], list[str]]]:
        entities: list[tuple[Literal["ligand", "interface"], list[str]]] = [
            ("ligand", list(self.cfg.flow.cluster_metrics))
        ]
        entities.append(("interface", ["interface_qcov"]))
        return entities

    @utils.ingest_flow_control
    def plan_clusters(self) -> None:
        from plinder.data.pipeline.score import plan_clustering

        for entity_type, metrics in self._cluster_entities():
            plan_clustering(
                self.plinder_dir,
                metrics=metrics,
                thresholds=list(self.cfg.flow.cluster_thresholds),
                source_batch_size=self.cfg.flow.symmetric_edge_source_batch_size,
                cover_batch_size=1,
                symmetric_bucket_count=self.cfg.flow.symmetric_edge_bucket_count,
                entity_type=entity_type,
            )

    @utils.ingest_flow_control
    def scatter_make_symmetric_edge_fragments(self) -> list[dict[str, Any]]:
        from plinder.data import clusters

        work: list[dict[str, Any]] = []
        for entity_type, _ in self._cluster_entities():
            plan = clusters.load_symmetric_edge_plan(
                self.plinder_dir,
                entity_type=entity_type,
            )
            work.extend(
                {"entity_type": entity_type, "batches": [batch]}
                for batch in plan["batches"]
            )
        return work or [{}]

    @utils.ingest_flow_control
    def make_symmetric_edge_fragments(self, work: dict[str, Any]) -> None:
        if not work:
            return
        tasks.make_symmetric_edge_fragments(
            data_dir=self.plinder_dir,
            batches=work["batches"],
            scratch_dir=Path(tempfile.gettempdir()) / "plinder-symmetric-fragments",
            threads=self.cfg.flow.clustering_cpu,
            force_update=self.cfg.data.force_update,
            entity_type=work["entity_type"],
        )

    @utils.ingest_flow_control
    def scatter_make_symmetric_edge_shards(self) -> list[dict[str, Any]]:
        from plinder.data import clusters

        work: list[dict[str, Any]] = []
        for entity_type, _ in self._cluster_entities():
            plan = clusters.load_symmetric_edge_plan(
                self.plinder_dir,
                entity_type=entity_type,
            )
            for metric in plan["metrics"]:
                for bucket in range(int(plan["bucket_count"])):
                    work.append(
                        {
                            "entity_type": entity_type,
                            "metric_buckets": [(str(metric), bucket)],
                        }
                    )
        return work or [{}]

    @utils.ingest_flow_control
    def make_symmetric_edge_shards(self, work: dict[str, Any]) -> None:
        if not work:
            return
        tasks.make_symmetric_edge_shards(
            data_dir=self.plinder_dir,
            metric_buckets=work["metric_buckets"],
            scratch_dir=Path(tempfile.gettempdir()) / "plinder-symmetric-edges",
            threads=self.cfg.flow.clustering_cpu,
            force_update=self.cfg.data.force_update,
            entity_type=work["entity_type"],
        )

    @utils.ingest_flow_control
    def scatter_make_component_reductions(self) -> list[Any]:
        work: list[dict[str, Any]] = []
        for entity_type, metrics in self._cluster_entities():
            batches = tasks.scatter_component_reduction_sources(
                data_dir=self.plinder_dir,
                metrics=metrics,
                batch_size=self.cfg.flow.component_reduction_source_batch_size,
                entity_type=entity_type,
            )
            work.extend(
                {"entity_type": entity_type, "metrics": metrics, "sources": batch}
                for batch in batches
                if batch
            )
        return work or [{}]

    @utils.ingest_flow_control
    def make_component_reductions(self, work: Any) -> None:
        if not work:
            return
        entity_type = work["entity_type"]
        metrics = work["metrics"]
        source_paths = work["sources"]
        tasks.make_component_reductions(
            data_dir=self.plinder_dir,
            source_paths=source_paths,
            metrics=metrics,
            thresholds=self.cfg.flow.cluster_thresholds,
            scratch_dir=Path(tempfile.gettempdir()) / "plinder-component-reductions",
            force_update=self.cfg.data.force_update,
            metric_workers=self.cfg.flow.component_reduction_metric_workers,
            entity_type=entity_type,
        )

    @utils.ingest_flow_control
    def merge_component_reductions(self) -> None:
        for entity_type, metrics in self._cluster_entities():
            tasks.merge_component_reductions(
                data_dir=self.plinder_dir,
                metrics=metrics,
                thresholds=self.cfg.flow.cluster_thresholds,
                entity_type=entity_type,
            )

    @utils.ingest_flow_control
    def scatter_make_set_covers(self) -> list[Any]:
        force_update = (
            self.cfg.data.force_update or self.cfg.flow.make_components_force_update
        )
        work: list[dict[str, Any]] = []
        for entity_type, metrics in self._cluster_entities():
            batches = tasks.scatter_make_set_covers(
                data_dir=self.plinder_dir,
                metrics=metrics,
                thresholds=self.cfg.flow.cluster_thresholds,
                stop_on_cluster=self.cfg.flow.make_components_stop_on_cluster,
                skip_existing_clusters=not force_update,
                entity_type=entity_type,
            )
            work.extend(
                {"entity_type": entity_type, "metric_threshold": batch}
                for batch in batches
                if batch
            )
        return work or [{}]

    @utils.ingest_flow_control
    def make_set_covers(self, work: Any) -> None:
        force_update = (
            self.cfg.data.force_update or self.cfg.flow.make_components_force_update
        )
        if not work:
            return
        entity_type = work["entity_type"]
        metric_thresholds = work["metric_threshold"]
        tasks.make_set_covers(
            data_dir=self.plinder_dir,
            metric_threshold=metric_thresholds,
            skip_existing_clusters=not force_update,
            scratch_dir=Path(tempfile.gettempdir()) / "plinder-set-covers",
            threads=self.cfg.flow.clustering_cpu,
            entity_type=entity_type,
        )

    @utils.ingest_flow_control
    def scatter_make_directed_set_covers(self) -> list[Any]:
        force_update = (
            self.cfg.data.force_update or self.cfg.flow.make_components_force_update
        )
        work: list[dict[str, Any]] = []
        for entity_type, metrics in self._cluster_entities():
            batches = tasks.scatter_make_directed_set_covers(
                data_dir=self.plinder_dir,
                metrics=metrics,
                thresholds=self.cfg.flow.cluster_thresholds,
                stop_on_cluster=self.cfg.flow.make_components_stop_on_cluster,
                skip_existing=not force_update,
                entity_type=entity_type,
            )
            work.extend(
                {"entity_type": entity_type, "metric_threshold": batch}
                for batch in batches
                if batch
            )
        return work or [{}]

    @utils.ingest_flow_control
    def make_directed_set_covers(self, work: Any) -> None:
        force_update = (
            self.cfg.data.force_update or self.cfg.flow.make_components_force_update
        )
        if not work:
            return
        entity_type = work["entity_type"]
        metric_thresholds = work["metric_threshold"]
        tasks.make_directed_set_covers(
            data_dir=self.plinder_dir,
            metric_threshold=metric_thresholds,
            skip_existing=not force_update,
            scratch_dir=Path(tempfile.gettempdir()) / "plinder-directed-set-covers",
            threads=self.cfg.flow.clustering_cpu,
            entity_type=entity_type,
        )

    @utils.ingest_flow_control
    def summarize_clusters(self) -> None:
        for entity_type, metrics in self._cluster_entities():
            tasks.summarize_clusters(
                data_dir=self.plinder_dir,
                metrics=metrics,
                thresholds=self.cfg.flow.cluster_thresholds,
                entity_type=entity_type,
            )

    @utils.ingest_flow_control
    def finalize_index(self) -> None:
        tasks.finalize_index(data_dir=self.plinder_dir)

    def run_stage(self, stage: str) -> None:
        """
        A stage is defined minimally as a {method} that
        is defined on the pipeline, optionally with
        a scatter_{method} to chunk inputs to feed into
        {method} and optionally a join_{method} to reduce
        the outputs produced by it.

        Parameters
        ----------
        stage : str
            name of the stage to run
        """
        scatter = getattr(self, f"scatter_{stage}", None)
        compute = getattr(self, stage)
        join = getattr(self, f"join_{stage}", None)
        chunks = None
        if scatter is not None:
            chunks = scatter()
        if chunks is not None:
            outs = [compute(chunk) for chunk in chunks]
        else:
            outs = [compute()]
        if join is not None:
            join(outs)

    def run(self) -> None:
        """
        Note that the order of operations defined here imply
        a directed acyclic graph. However, we can manually add
        a cycle to the graph by feeding the output of
        join_make_entries back into make_entries, e.g. in order
        to bump up resource requests for systems that failed to
        run successfully on lower memory pods. The actual step
        execution is offloaded to a workflow orchestration
        framework in this case, and the run method is not called
        directly. However, this is an exceedingly convenient way
        to run the pipeline locally.
        """
        for stage in tasks.STAGES:
            self.run_stage(stage)


if __name__ == "__main__":
    pipe = IngestPipeline()
    pipe.run()
