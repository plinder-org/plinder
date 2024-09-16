# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
from pathlib import Path
from typing import Any, Optional

from omegaconf import DictConfig

from plinder.core.utils.log import setup_logger
from plinder.data.pipeline import config, tasks, utils

LOG = setup_logger(__name__)


class IngestPipeline:
    """
    Mimic the required metaflow DAG pattern of
        - scatter
        - compute
        - join
    by convention in method names. This could be
    enforced with a metaclass construct but not enough time.

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
    def make_dbs(self) -> None:
        tasks.make_dbs(
            data_dir=self.plinder_dir,
            sub_databases=self.cfg.scorer.sub_databases,
        )

    @utils.ingest_flow_control
    def scatter_make_entries(self) -> list[list[str]]:
        force_update = (
            self.cfg.data.force_update or self.cfg.flow.make_entries_force_update
        )
        chunks: list[list[str]] = tasks.scatter_make_entries(
            data_dir=self.plinder_dir,
            batch_size=self.cfg.flow.make_entries_batch_size,
            two_char_codes=self.cfg.context.two_char_codes,
            pdb_ids=self.cfg.context.pdb_ids,
            force_update=force_update,
        )
        return chunks

    @utils.ingest_flow_control
    def make_entries(self, pdb_dirs: list[str]) -> list[str]:
        force_update = (
            self.cfg.data.force_update or self.cfg.flow.make_entries_force_update
        )
        failed: list[str] = tasks.make_entries(
            data_dir=self.plinder_dir,
            pdb_dirs=pdb_dirs,
            force_update=force_update,
            cpu=self.cfg.flow.make_entries_cpu,
            annotation_cfg=self.cfg.annotation,
            entry_cfg=self.cfg.entry,
        )
        return failed

    @utils.ingest_flow_control
    def join_make_entries(self, reruns: list[list[str]]) -> list[str]:
        catted = []
        for rerun in reruns:
            catted.extend([item[-4:] for item in rerun])
        return catted

    @utils.ingest_flow_control
    def scatter_structure_qc(self) -> list[list[str]]:
        chunks: list[list[str]] = tasks.scatter_structure_qc(
            data_dir=self.plinder_dir,
            two_char_codes=self.cfg.context.two_char_codes,
            batch_size=self.cfg.flow.download_rcsb_files_batch_size,
        )
        return chunks

    @utils.ingest_flow_control
    def structure_qc(self, two_char_codes: list[str]) -> None:
        tasks.structure_qc(
            data_dir=self.plinder_dir,
            two_char_codes=two_char_codes,
        )

    @utils.ingest_flow_control
    def join_structure_qc(self, qcs: list[None]) -> None:
        force_update = (
            self.cfg.data.force_update or self.cfg.flow.structure_qc_force_update
        )
        utils.create_index(
            data_dir=self.plinder_dir,
            force_update=force_update,
        )

    @utils.ingest_flow_control
    def scatter_make_ligands(self) -> list[list[str]]:
        chunks: list[list[str]] = tasks.scatter_make_ligands(
            data_dir=self.plinder_dir,
            batch_size=self.cfg.flow.make_ligands_batch_size,
            two_char_codes=self.cfg.context.two_char_codes,
            pdb_ids=self.cfg.context.pdb_ids,
        )
        return chunks

    @utils.ingest_flow_control
    def make_ligands(self, pdb_ids: list[str]) -> None:
        tasks.make_ligands(data_dir=self.plinder_dir, pdb_ids=pdb_ids)

    @utils.ingest_flow_control
    def compute_ligand_fingerprints(self) -> None:
        tasks.compute_ligand_fingerprints(
            data_dir=self.plinder_dir,
            split_char=self.cfg.ligand.ligand_id_split_char,
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
            save_top_k_similar_ligands=self.cfg.ligand.save_top_k_similar_ligands,
            multiply_by=self.cfg.ligand.multiply_by,
            number_id_col=self.cfg.ligand.number_id_col,
        )

    @utils.ingest_flow_control
    def make_mmp_index(self) -> None:
        tasks.make_mmp_index(data_dir=self.plinder_dir)

    @utils.ingest_flow_control
    def scatter_make_system_archives(self) -> list[list[str]]:
        chunks: list[list[str]] = tasks.scatter_make_system_archives(
            data_dir=self.plinder_dir,
            two_char_codes=self.cfg.context.two_char_codes,
            batch_size=self.cfg.flow.download_rcsb_files_batch_size,
        )
        return chunks

    @utils.ingest_flow_control
    def make_system_archives(self, two_char_codes: list[str]) -> None:
        tasks.make_system_archives(
            data_dir=self.plinder_dir,
            two_char_codes=two_char_codes,
        )

    @utils.ingest_flow_control
    def make_sub_dbs(self) -> None:
        tasks.make_sub_dbs(
            data_dir=self.plinder_dir,
            sub_databases=self.cfg.scorer.sub_databases,
        )

    @utils.ingest_flow_control
    def scatter_run_batch_searches(self) -> list[list[str]]:
        chunks: list[list[str]] = tasks.scatter_protein_scoring(
            data_dir=self.plinder_dir,
            batch_size=self.cfg.flow.run_batch_searches_batch_size,
            two_char_codes=self.cfg.context.two_char_codes,
            pdb_ids=self.cfg.context.pdb_ids,
        )
        return chunks

    @utils.ingest_flow_control
    def run_batch_searches(self, pdb_ids: list[str]) -> None:
        tasks.run_batch_searches(
            data_dir=self.plinder_dir,
            pdb_ids=pdb_ids,
            scorer_cfg=self.cfg.scorer,
        )

    @utils.ingest_flow_control
    def scatter_make_batch_scores(self) -> list[list[str]]:
        chunks: list[list[str]] = tasks.scatter_missing_scores(
            data_dir=self.plinder_dir,
            batch_size=self.cfg.flow.make_batch_scores_batch_size,
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
        )

    @utils.ingest_flow_control
    def scatter_collate_partitions(self) -> list[list[str]]:
        chunks: list[list[str]] = tasks.scatter_collate_partitions()
        return chunks

    @utils.ingest_flow_control
    def collate_partitions(self, partition: list[str]) -> None:
        tasks.collate_partitions(data_dir=self.plinder_dir, partition=partition)

    @utils.ingest_flow_control
    def scatter_make_components_and_communities(self) -> list[list[tuple[str, int]]]:
        force_update = (
            self.cfg.data.force_update or self.cfg.flow.make_components_force_update
        )
        chunks: list[
            list[tuple[str, int]]
        ] = tasks.scatter_make_components_and_communities(
            data_dir=self.plinder_dir,
            metrics=self.cfg.flow.cluster_metrics,
            thresholds=self.cfg.flow.cluster_thresholds,
            stop_on_cluster=self.cfg.flow.make_components_stop_on_cluster,
            skip_existing_clusters=not force_update,
        )
        return chunks

    @utils.ingest_flow_control
    def make_components_and_communities(
        self, metric_thresholds: list[tuple[str, int]]
    ) -> None:
        force_update = (
            self.cfg.data.force_update or self.cfg.flow.make_components_force_update
        )
        tasks.make_components_and_communities(
            data_dir=self.plinder_dir,
            metric_threshold=metric_thresholds,
            skip_existing_clusters=not force_update,
        )

    @utils.ingest_flow_control
    def scatter_make_splits(self) -> list[list[tuple[DictConfig, str]]]:
        chunks: list[list[tuple[DictConfig, str]]] = tasks.scatter_make_splits(
            data_dir=self.plinder_dir,
            split_config_dir=self.cfg.flow.split_config_dir,
        )
        return chunks

    @utils.ingest_flow_control
    def make_splits(self, cfg_and_path: list[tuple[DictConfig, str]]) -> None:
        tasks.make_splits(data_dir=self.plinder_dir, cfg_and_path=cfg_and_path)

    @utils.ingest_flow_control
    def scatter_compute_ligand_leakage(self) -> list[list[tuple[str, str, str]]]:
        chunks: list[list[tuple[str, str, str]]] = tasks.scatter_compute_ligand_leakage(
            data_dir=self.plinder_dir,
            test_leakage=self.cfg.flow.test_leakage,
        )
        return chunks

    @utils.ingest_flow_control
    def compute_ligand_leakage(self, inputs: list[tuple[str, str, str]]) -> None:
        tasks.compute_ligand_leakage(
            data_dir=self.plinder_dir,
            inputs=inputs,
        )

    @utils.ingest_flow_control
    def scatter_compute_protein_leakage(self) -> list[list[tuple[str, str, str]]]:
        chunks: list[
            list[tuple[str, str, str]]
        ] = tasks.scatter_compute_protein_leakage(
            data_dir=self.plinder_dir,
            test_leakage=self.cfg.flow.test_leakage,
        )
        return chunks

    @utils.ingest_flow_control
    def compute_protein_leakage(self, inputs: list[tuple[str, str, str]]) -> None:
        tasks.compute_protein_leakage(
            data_dir=self.plinder_dir,
            inputs=inputs,
        )

    @utils.ingest_flow_control
    def scatter_make_links(self) -> list[list[str]]:
        chunks: list[list[str]] = tasks.scatter_make_links(
            data_dir=self.plinder_dir,
            search_dbs=self.cfg.scorer.sub_databases,
        )
        return chunks

    @utils.ingest_flow_control
    def make_links(self, search_dbs: list[str]) -> None:
        tasks.make_links(
            data_dir=self.plinder_dir,
            search_dbs=search_dbs,
        )

    @utils.ingest_flow_control
    def make_linked_structures(self) -> None:
        force_update = (
            self.cfg.data.force_update
            or self.cfg.flow.make_linked_structures_force_update
        )
        tasks.make_linked_structures(
            data_dir=self.plinder_dir,
            search_dbs=self.cfg.flow.sub_databases,
            cpu=self.cfg.flow.make_linked_structures_cpu,
            force_update=force_update,
        )

    @utils.ingest_flow_control
    def scatter_score_linked_structures(self) -> list[list[tuple[str, str]]]:
        chunks: list[list[tuple[str, str]]] = tasks.scatter_score_linked_structures(
            data_dir=self.plinder_dir,
            search_dbs=self.cfg.flow.sub_databases,
            batch_size=self.cfg.flow.score_linked_structures_batch_size,
        )
        return chunks

    @utils.ingest_flow_control
    def score_linked_structures(self, system_ids: list[tuple[str, str]]) -> None:
        force_update = (
            self.cfg.data.force_update
            or self.cfg.flow.score_linked_structures_force_update
        )
        tasks.score_linked_structures(
            data_dir=self.plinder_dir,
            search_dbs=self.cfg.flow.sub_databases,
            system_ids=system_ids,
            cpu=self.cfg.flow.score_linked_structures_cpu,
            force_update=force_update,
        )
        return

    @utils.ingest_flow_control
    def join_score_linked_structures(self, outputs: list[None]) -> None:
        utils.mp_pack_linked_structures(data_dir=self.plinder_dir, structures=False)
        utils.consolidate_linked_scores(data_dir=self.plinder_dir)

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
