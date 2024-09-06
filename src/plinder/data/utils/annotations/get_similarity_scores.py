# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
from __future__ import annotations

import shutil
import subprocess
from collections import abc, defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

import biotite.sequence as seq
import biotite.sequence.align as align
import numpy as np
import pandas as pd
import pyarrow
import pyarrow.parquet as pq
from pyarrow import csv
from tqdm import tqdm

from plinder.core.utils.log import setup_logger
from plinder.data import databases
from plinder.data.pipeline.config import FoldseekConfig, MMSeqsConfig
from plinder.data.pipeline.utils import load_entries_from_zips
from plinder.data.utils.annotations.aggregate_annotations import Entry, System

LOG = setup_logger(__name__)

SORT_ORDER = [
    ("similarity", "descending"),
    ("query_system", "ascending"),
    ("target_system", "ascending"),
]
INFO_COLUMNS = (
    "query_system",
    "target_system",
    "protein_mapping",
    "protein_mapper",
)
SCORE_NAMES = (
    "protein_lddt",
    "protein_lddt_qcov",
    "protein_qcov",
    "protein_fident",
    "protein_fident_qcov",
    "protein_seqsim",
    "protein_seqsim_qcov",
    "pocket_lddt",
    "pocket_lddt_qcov",
    "pocket_qcov",
    "pocket_fident",
    "pocket_fident_qcov",
    "pli_qcov",
    "pli_unique_qcov",
)

_ChainInstanceMapping = str
_ChainPairType = tuple[_ChainInstanceMapping, _ChainInstanceMapping]
_SimilarityScoreDictType = dict[str, float]


def get_sequence_similarity(seq_str1: str, seq_str2: str) -> tuple[float, float]:
    """
    Calculate the similarity score of an alignment.

    If the alignment contains more than two sequences,
    all pairwise scores are counted.

    Parameters
    ----------
    seq_str1 : str
        First protein sequence
    seq_str2 : str
        Second protein sequence

    Returns
    -------
    tuple[float, float]
        Sequence identity and sequence similarity score.
    """
    non_canonical_aa: dict[str, str | int | None] = {
        "X": "A",  # Replace X (any a.a with alanine)
        "B": "D",  # Replace Asx ( with aspartic acid)
        "J": "L",  # Replace Xle ( with leucine)
        "Z": "E",  # Replace Glx (any a.a with glutamic acid)
        "U": "C",  # Replace Selenocysteine(sec) (with cysteine)
        "O": "K",  # replace Pyrrolysine(Pyl) (any a.a with alanine)
    }
    seq_str1 = seq_str1.translate(str.maketrans(non_canonical_aa))
    seq_str2 = seq_str2.translate(str.maketrans(non_canonical_aa))
    seq1_arr = np.array(list(seq_str1))
    seq2_arr = np.array(list(seq_str2))
    gap1 = np.where(seq1_arr == "-")
    gap2 = np.where(seq2_arr == "-")
    all_gaps = np.concatenate([gap1[0], gap2[0]])
    seq1_arr = np.delete(seq1_arr, all_gaps)
    seq2_arr = np.delete(seq2_arr, all_gaps)
    ungapped_seq_str1 = "".join(list(seq1_arr))
    ungapped_seq_str2 = "".join(list(seq2_arr))
    seq1 = seq.ProteinSequence(ungapped_seq_str1)
    seq2 = seq.ProteinSequence(ungapped_seq_str2)
    matrix = align.SubstitutionMatrix.std_protein_matrix()
    trace = align.Alignment.trace_from_strings([ungapped_seq_str1, ungapped_seq_str2])
    ali = align.Alignment([seq1, seq2], trace)
    codes = align.alignment.get_codes(ali)
    matrix = matrix.score_matrix()

    # Sum similarity scores (without gaps)
    scores = []

    # Iterate over all positions
    for pos in range(codes.shape[1]):
        column = codes[:, pos]
        # Iterate over all possible pairs
        # Do not count self-similarity
        # and do not count similarity twice (not S(i,j) and S(j,i))
        for i in range(codes.shape[0]):
            for j in range(i + 1, codes.shape[0]):
                code_i = column[i]
                code_j = column[j]
                # Ignore gaps
                if code_i != -1 and code_j != -1:
                    tmp_score = max(0, matrix[code_i, code_j])
                    normalizing_const = max(
                        matrix[code_i, code_i], matrix[code_j, code_j]
                    )
                    scores.append(tmp_score / normalizing_const > 0.2)
    similarity_score = sum(scores) / len(scores)
    return (align.get_sequence_identity(ali), similarity_score)


def get_sequence_similarity_helper(seq_str1: str, seq_str2: str) -> float:
    try:
        return get_sequence_similarity(seq_str1, seq_str2)[1]
    except Exception:
        return 0


def run_alignment(
    aln_type: str,
    query_db: Path,
    target_db: Path,
    search_db: Path,
    aln_file: Path,
    alignment_config: FoldseekConfig | MMSeqsConfig,
    tmp_dir: Path = Path.cwd() / "tmp",
    remove_tmp: bool = True,
) -> None:
    if search_db.with_suffix(".dbtype").exists():
        search_db.with_suffix(".dbtype").unlink()

    search_commands = [
        aln_type,
        "search",
        str(query_db),
        str(target_db),
        str(search_db),
        str(tmp_dir),
        "-a",
        "-e",
        f"{alignment_config.evalue}",
        "-s",
        f"{alignment_config.sensitivity}",
        "--max-seqs",
        f"{alignment_config.max_seqs}",
        "-c",
        f"{alignment_config.coverage}",
        "--cov-mode",
        "2",
        "--min-seq-id",
        f"{alignment_config.min_seq_id}",
    ]
    if aln_type == "foldseek":
        search_commands += ["--sort-by-structure-bits", "0"]
    convert_commands = [
        aln_type,
        "convertalis",
        str(query_db),
        str(target_db),
        str(search_db),
        str(aln_file.with_suffix(".tsv")),
        "--format-mode",
        "4",
        "--format-output",
        "query,target,qlen,fident,alnlen,qstart,qend,tstart,tend,evalue,bits,qcov,tcov,qaln,taln",
    ]
    if aln_type == "foldseek":
        convert_commands[-1] += ",lddt,lddtfull"
    subprocess.check_call(search_commands, stdout=subprocess.DEVNULL)
    subprocess.check_call(convert_commands, stdout=subprocess.DEVNULL)

    table = csv.read_csv(
        aln_file.with_suffix(".tsv"), parse_options=csv.ParseOptions(delimiter="\t")
    )
    table = table.append_column(
        "query_pdb_id",
        pyarrow.array(
            [
                str(x)
                .replace("_xyz-enrich.cif.gz", "")
                .replace("_xyz-enrich.cif", "")
                .replace("pdb_0000", "")[:4]
                for x in table["query"]
            ]
        ),
    )
    if search_db != "pred":
        table = table.append_column(
            "target_pdb_id",
            pyarrow.array(
                [
                    str(x)
                    .replace("_xyz-enrich.cif.gz", "")
                    .replace("_xyz-enrich.cif", "")
                    .replace("pdb_0000", "")[:4]
                    for x in table["target"]
                ]
            ),
        )
    if aln_file.with_suffix(".parquet").exists():
        shutil.rmtree(aln_file.with_suffix(".parquet"))
    pq.write_to_dataset(
        table, aln_file.with_suffix(".parquet"), partition_cols=["query_pdb_id"]
    )

    # Cleanup
    if remove_tmp and tmp_dir.exists():
        shutil.rmtree(tmp_dir)
        aln_file.with_suffix(".tsv").unlink()
    for filename in search_db.parent.glob(f"{search_db.name}*"):
        if filename.is_dir():
            shutil.rmtree(filename)
        elif filename.is_file():
            filename.unlink()


def combine_scores(
    q_t_scores: _SimilarityScoreDictType,
    q_t_mappings: dict[str, list[_ChainPairType]],
    protein_chain_mapper: str = "",
) -> dict[str, str | float]:
    def make_mapping(mapping: list[_ChainPairType]) -> str:
        return ";".join(f"{a}:{b}" if len(b) else a for a, b in mapping)

    q_t_scores_combined: dict[str, str | float] = {}
    for s in SCORE_NAMES:
        for suffix in ["_max", "_weighted_max", "_weighted_sum", ""]:
            if (
                f"{s}_foldseek{suffix}" not in q_t_scores
                and f"{s}_mmseqs{suffix}" not in q_t_scores
            ):
                continue
            s1, s2 = (
                float(q_t_scores.get(f"{s}_foldseek{suffix}", 0)),
                float(q_t_scores.get(f"{s}_mmseqs{suffix}", 0)),
            )
            source = "foldseek"
            if np.abs(s1 - s2) < 1e-6:
                q_t_scores_combined[f"{s}{suffix}"] = s1
                q_t_scores_combined[f"{s}{suffix}_source"] = "both"
            else:
                index = np.argmax([s1, s2])
                q_t_scores_combined[f"{s}{suffix}"] = [s1, s2][index]
                source = ["foldseek", "mmseqs"][index]
                q_t_scores_combined[f"{s}{suffix}_source"] = source
            source_name = "mmseqs" if source == "mmseqs" else "foldseek"
            n = f"{s}_{source_name}{suffix}"
            if n in q_t_mappings:
                q_t_scores_combined[f"{s}{suffix}_mapping"] = make_mapping(
                    q_t_mappings[n]
                )
    if protein_chain_mapper != "":
        q_t_scores_combined["protein_mapping"] = make_mapping(
            q_t_mappings.get(f"{protein_chain_mapper}_weighted_sum", [])
        )
        q_t_scores_combined["protein_mapper"] = (
            "foldseek" if "lddt" in protein_chain_mapper else "mmseqs"
        )
    return q_t_scores_combined


@dataclass
class Scorer:
    entries: dict[str, Entry]
    source_to_full_db_file: dict[str, Path]
    db_dir: Path
    scores_dir: Path
    protein_chain_mappers: list[str] = field(
        default_factory=lambda: [
            "protein_lddt_qcov_foldseek",
            "protein_fident_qcov_mmseqs",
        ]
    )
    ligand_chain_mappers: list[str] = field(
        default_factory=lambda: [
            "pocket_lddt_qcov_foldseek",
            "pocket_fident_qcov_mmseqs",
        ]
    )
    # default minimum threshold for every metric BEFORE scaling to 100
    minimum_threshold: int = 0
    # optional custom minimum threshold for each metric BEFORE scaling to 100
    minimum_thresholds: dict[str, int] = field(
        default_factory=lambda: {"pli_qcov": 0, "pocket_qcov": 0, "pli_unique_qcov": 0}
    )

    def __post_init__(self) -> None:
        self.db_dir.mkdir(exist_ok=True, parents=True)
        self.scores_dir.mkdir(exist_ok=True, parents=True)

    def make_dbs(self) -> None:
        databases.make_sub_dbs(self.db_dir, self.source_to_full_db_file, self.entries)

    @staticmethod
    def get_config(search_db: str, aln_type: str) -> FoldseekConfig | MMSeqsConfig:
        config: FoldseekConfig | MMSeqsConfig
        if aln_type == "foldseek":
            config = FoldseekConfig()
        else:
            config = MMSeqsConfig()
        if search_db in ["apo", "pred"]:
            config.coverage = 0.9
            config.min_seq_id = 0.9
        return config

    def run_alignments(
        self,
        entry_ids: list[str],
        search_db: str,
        output_folder: Path,
        overwrite: bool = False,
    ) -> None:
        output_folder.mkdir(exist_ok=True)
        for aln_type in ["mmseqs", "foldseek"]:
            sub_db = output_folder / search_db / aln_type
            sub_db.mkdir(exist_ok=True, parents=True)
            db_ids = databases.get_db_ids(
                self.entries, "holo", aln_type, entry_ids=entry_ids
            )
            databases.make_sub_db(
                db_ids,
                self.source_to_full_db_file[f"holo_{aln_type}"],
                sub_db,
                aln_type,
            )
            tmp_dir = sub_db / f"tmp_{search_db}_{aln_type}"
            tmp_dir.mkdir(exist_ok=True, parents=True)
            aln_file = sub_db / f"aln_{search_db}.tsv"
            LOG.info("run_alignments calling run_alignment:")
            LOG.info(f"    sub_db={sub_db}")
            LOG.info(f"    aln_file={aln_file.with_suffix('.tsv')}")
            LOG.info(f"    aln_type={aln_type}")
            try:
                run_alignment(
                    aln_type=aln_type,
                    query_db=sub_db / sub_db.name,
                    target_db=self.db_dir
                    / f"{search_db}_{aln_type}"
                    / f"{search_db}_{aln_type}",
                    search_db=sub_db / "search",
                    aln_file=aln_file.with_suffix(".tsv"),
                    tmp_dir=tmp_dir / output_folder.stem,
                    alignment_config=self.get_config(search_db, aln_type),
                )
            except Exception as e:
                scratch = (
                    Path(*output_folder.parts[:3])
                    / "scratch"
                    / "scores"
                    / "run_alignment_failures"
                )
                scratch.mkdir(exist_ok=True, parents=True)
                (scratch / f"{search_db}_{aln_type}.txt").write_text(f"{repr(e)}: {e}")
                LOG.error(f"scoring: Error for {search_db}_{aln_type}: {e}")
                continue
            aln_dir = self.db_dir / f"{search_db}_{aln_type}" / "aln"
            aln_dir.mkdir(exist_ok=True, parents=True)
            for pdb_id in tqdm(entry_ids):
                pdb_id_file = (
                    aln_file.with_suffix(".parquet") / f"query_pdb_id={pdb_id}"
                )
                if not pdb_id_file.exists():
                    scratch = (
                        Path(*output_folder.parts[:3])
                        / "scratch"
                        / "scores"
                        / "run_alignments_failures"
                    )
                    scratch.mkdir(exist_ok=True, parents=True)
                    (scratch / f"{search_db}_{aln_type}_{pdb_id}.txt").write_text("")
                    continue
                pdb_id_df = pd.read_parquet(pdb_id_file)
                if not pdb_id_df.empty:
                    pdb_id_df.to_parquet(aln_dir / f"{pdb_id}.parquet")
                else:
                    scratch = (
                        Path(*output_folder.parts[:3])
                        / "scratch"
                        / "scores"
                        / "run_alignments_empty"
                    )
                    scratch.mkdir(exist_ok=True, parents=True)
                    (scratch / f"{search_db}_{aln_type}_{pdb_id}.txt").write_text("")

    def get_score_df(
        self, data_dir: Path, pdb_id: str, search_db: str, overwrite: bool = True
    ) -> Path:
        """
        Convert aligmnent results to mapped alignment results. Then
        aggregate the mapped alignment results into the scores dataset.
        """
        score_df_path = self.db_dir / f"search_db={search_db}" / f"{pdb_id}.parquet"
        if overwrite or not score_df_path.exists():
            LOG.info(f"get_score_df: aggregating scores for {pdb_id} to {search_db}")
        else:
            LOG.info(f"get_score_df: skipping existing {score_df_path}")
            return score_df_path
        for aln_type in ["foldseek", "mmseqs"]:
            pdb_id_file = (
                self.db_dir / f"{search_db}_{aln_type}" / "aln" / f"{pdb_id}.parquet"
            )
            if not pdb_id_file.exists():
                LOG.info(f"get_score_df: pdb_id_file={pdb_id_file} does not exist")
                continue
            pdb_file = (
                self.db_dir
                / f"{search_db}_{aln_type}"
                / "mapped_aln"
                / f"{pdb_id}.parquet"
            )
            pdb_file.parent.mkdir(exist_ok=True, parents=True)
            if overwrite or not pdb_file.exists():
                self.entries = {}
                entries_to_load = {pdb_id}
                if search_db != "pred" and pdb_id_file.exists():
                    entries_to_load |= set(
                        pd.read_parquet(pdb_id_file, columns=["target_pdb_id"])[
                            "target_pdb_id"
                        ]
                    )
                entries_to_load = entries_to_load.difference(self.entries.keys())
                LOG.info(f"entries_to_load pdb_id={pdb_id} {len(entries_to_load)}")
                LOG.info(
                    f"loading {len(entries_to_load)} (additional) entries for {pdb_id}"
                )
                self.entries.update(
                    load_entries_from_zips(
                        data_dir=data_dir,
                        pdb_ids=entries_to_load,
                        load_for_scoring=True,
                    )
                )

                try:
                    LOG.info(
                        f"mapping aligmnet df for {pdb_id} to {search_db} for {aln_type}"
                    )
                    self.map_alignment_df(pdb_id_file, aln_type, search_db).to_parquet(
                        pdb_file, index=True
                    )
                except Exception as e:
                    scratch = (
                        Path(*pdb_id_file.parts[:3])
                        / "scratch"
                        / "scores"
                        / "map_alignment_df_failures"
                    )
                    scratch.mkdir(exist_ok=True, parents=True)
                    (scratch / f"{search_db}_{aln_type}_{pdb_id}.txt").write_text(
                        f"{repr(e)}: {e}"
                    )
                    LOG.error(
                        f"scoring: Error in map_alignment_df: {pdb_id} searching against {search_db} with {aln_type}: {repr(e)}"
                    )
                    continue
            else:
                LOG.info(f"skipping creating {pdb_file} because it already exists")

        try:
            score_df_path.parent.mkdir(exist_ok=True, parents=True)
            LOG.info(f"aggregating scores for {pdb_id} to {search_db}")
            df = self.aggregate_scores(pdb_id, search_db=search_db)
            if df is not None and not df.empty:
                df.to_parquet(score_df_path, index=False)
        except Exception as e:
            scratch = (
                Path(*score_df_path.parts[:3])
                / "scratch"
                / "scores"
                / "aggregate_scores_failures"
            )
            scratch.mkdir(exist_ok=True, parents=True)
            (scratch / f"{search_db}_{pdb_id}.txt").write_text(f"{repr(e)}: {e}")
            LOG.error(
                f"scoring: Error in aggregate_scores: {pdb_id} searching against {search_db}: {repr(e)}"
            )
        return score_df_path

    def load_alignments(
        self,
        source_to_aln_file: dict[str, Path],
        search_db: str = "holo",
    ) -> pd.DataFrame:
        data = []
        for source, aln_file in source_to_aln_file.items():
            sdb, aln_type = source.split("_")
            if sdb != search_db:
                continue
            if not aln_file.exists():
                continue
            aln_df = pd.read_parquet(aln_file)
            aln_df["source"] = aln_type
            aln_df["qrnum"] = aln_df["qrnum"].apply(
                lambda x: dict([(int(i), int(r)) for i, r in x])
            )
            aln_df["trnum"] = aln_df["trnum"].apply(
                lambda x: dict([(int(i), int(r)) for i, r in x])
            )
            if "lddtfull" in aln_df.columns:
                aln_df["lddtfull"] = aln_df["lddtfull"].apply(
                    lambda x: dict([(int(i), l) for i, l in x])
                )
            data.append(aln_df)
        if len(data):
            df = pd.concat(data)
            return df.sort_index()
        else:
            return pd.DataFrame()

    def map_alignment_df(
        self, df_file: Path, aln_type: str, search_db: str
    ) -> pd.DataFrame:
        df = pd.read_parquet(df_file)
        if aln_type == "foldseek":
            df["query"] = df["query"].replace(
                {"_xyz-enrich.cif.gz": "", "_xyz-enrich.cif": "", "pdb_0000": ""},
                regex=True,
            )
            if search_db == "pred":
                df["target"] = df["target"].replace(
                    {"-F1-model_v4.cif": "", "-F1-model_v4": "", "AF-": ""}, regex=True
                )
            else:
                df["target"] = df["target"].replace(
                    {"_xyz-enrich.cif.gz": "", "_xyz-enrich.cif": "", "pdb_0000": ""},
                    regex=True,
                )
        df["query_chain_mapped"] = (
            df["query"]
            .str.split("_", expand=True)
            .apply(lambda x: self.entries[x[0]].author_to_asym.get(x[1], None), axis=1)
        )
        df = df.dropna(subset=["query_chain_mapped"]).reset_index(drop=True)
        df["query_entry"] = df["query"].str.split("_", expand=True)[0]
        if search_db == "pred":
            df["target_chain_mapped"] = "A"
        else:
            df["target_chain_mapped"] = (
                df["target"]
                .str.split("_", expand=True)
                .apply(
                    lambda x: self.entries[x[0]].author_to_asym.get(x[1], None)
                    if x[0] in self.entries
                    else None,
                    axis=1,
                )
            )
            df = df.dropna(subset=["target_chain_mapped"]).reset_index(drop=True)
        df["target_entry"] = df["target"].str.split("_", expand=True)[0]
        df["qaln"] = df["qaln"].str.upper()
        df["taln"] = df["taln"].str.upper()
        df["seqsim"] = df[["qaln", "taln"]].apply(
            lambda x: get_sequence_similarity_helper(x["qaln"], x["taln"]), axis=1
        )
        df["seqsim_qcov"] = df["seqsim"] * df["qcov"]
        df["fident_qcov"] = df["fident"] * df["qcov"]
        if aln_type == "foldseek":
            df["lddt_qcov"] = df["lddt"] * df["qcov"]
            df["lddtaln"] = df["lddtfull"].apply(lambda x: x.split(","))
        df = df.apply(
            lambda x: self.map_row(x, aln_type=aln_type, search_db=search_db), axis=1
        )
        if aln_type == "foldseek":
            df.drop(columns=["lddtaln"], inplace=True)
        df["source"] = aln_type
        df.set_index(
            [
                "query_entry",
                "target_entry",
                "query_chain_mapped",
                "target_chain_mapped",
                "source",
            ],
            inplace=True,
        )
        return df

    def map_row(self, parts: pd.Series, aln_type: str, search_db: str) -> pd.Series:
        parts["qrnum"] = []
        parts["trnum"] = []
        parts["lddtfull"] = []
        q_i, t_i = int(parts["qstart"]) - 1, int(parts["tstart"]) - 1
        aln_index = 0
        for x, (q_a, t_a) in enumerate(zip(parts["qaln"], parts["taln"])):
            if q_a != "-" and t_a != "-":
                q_chain = self.entries[parts["query_entry"]].chains[
                    parts["query_chain_mapped"]
                ]
                t_chain = (
                    self.entries[parts["target_entry"]].chains[
                        parts["target_chain_mapped"]
                    ]
                    if search_db != "pred"
                    else None
                )
                if aln_type == "foldseek":
                    # it's the residue_index, map to residue_number
                    q_n = q_chain.residue_index_to_number.get(q_i, None)
                    if t_chain is not None:
                        t_n = t_chain.residue_index_to_number.get(t_i, None)
                    else:
                        t_n = None
                    if q_n is not None:
                        parts["qrnum"].append((x, q_n))
                        parts["lddtfull"].append(
                            (x, float(parts["lddtaln"][aln_index]))
                        )
                    if t_n is not None:
                        parts["trnum"].append((x, t_n))
                elif aln_type == "mmseqs":
                    if q_i in q_chain.residues:
                        parts["qrnum"].append((x, q_i))
                    if t_chain is not None and t_i in t_chain.residues:
                        parts["trnum"].append((x, t_i))
                aln_index += 1
            if q_a != "-":
                q_i += 1
            if t_a != "-":
                t_i += 1
        return parts

    def get_protein_scores_pair(self, aln: pd.DataFrame) -> dict[str, float]:
        scores = {}
        for source, data in aln.iterrows():
            for score in ["qcov", "fident", "seqsim", "fident_qcov", "seqsim_qcov"]:
                scores[f"protein_{score}_{source}"] = data[score]
            if source == "foldseek":
                for score in ["lddt", "lddt_qcov"]:
                    scores[f"protein_{score}_{source}"] = data[score]
        return scores

    def get_protein_scores(
        self,
        query_target_entry_alignments: pd.DataFrame,
        query_system: System,
        target_protein_chains: list[str],
        query_system_length: int,
    ) -> tuple[
        dict[str, list[_ChainPairType]],
        _SimilarityScoreDictType,
        dict[_ChainPairType, pd.DataFrame],
        str,
    ]:
        scores: _SimilarityScoreDictType = defaultdict(float)
        mappings: dict[str, list[_ChainPairType]] = defaultdict(list)
        max_chain_lengths: dict[str, float] = defaultdict(float)
        protein_chain_mapper = ""
        s_matrix = np.zeros(
            (len(query_system.protein_chains_asym_id), len(target_protein_chains))
        )
        for i, q_instance_chain in enumerate(query_system.protein_chains_asym_id):
            q_chain = q_instance_chain.split(".")[1]
            q_chain_length = self.entries[query_system.pdb_id].chains[q_chain].length
            for j, t_instance_chain in enumerate(target_protein_chains):
                t_chain = t_instance_chain.split(".")[1]
                try:
                    aln = query_target_entry_alignments.loc[q_chain].loc[t_chain]
                except KeyError:
                    LOG.debug("get_protein_scores indexing failed, aln is None")
                    aln = None
                if aln is not None:
                    pair = (q_instance_chain, t_instance_chain)
                    pair_scores = self.get_protein_scores_pair(aln)
                    if protein_chain_mapper == "":
                        for chain_mapper in self.protein_chain_mappers:
                            if chain_mapper in pair_scores:
                                protein_chain_mapper = chain_mapper
                                break
                    s_matrix[i, j] = pair_scores.get(protein_chain_mapper, 0)
                    for score in pair_scores:
                        if (
                            pair_scores[score] * q_chain_length
                            > scores[f"{score}_weighted_max"]
                        ):
                            scores[f"{score}_weighted_max"] = (
                                pair_scores[score] * q_chain_length
                            )
                            max_chain_lengths[f"{score}_weighted_max"] = q_chain_length
                            mappings[f"{score}_weighted_max"] = [pair]
                        if pair_scores[score] > scores[f"{score}_max"]:
                            scores[f"{score}_max"] = pair_scores[score]
                            mappings[f"{score}_max"] = [pair]
        for score in scores:
            if score.endswith("_weighted_max") and max_chain_lengths[score] > 0:
                scores[score] /= max_chain_lengths[score]

        # do chain mapping
        # Find the highest scores along with their indices in the matrix
        alns = {}
        while np.any(s_matrix):
            # changed score -> score_tag conflicting variable name
            score_val = np.amax(s_matrix)
            if score_val == 0:
                break
            q_idx, t_idx = np.unravel_index(np.argmax(s_matrix), s_matrix.shape)
            q_instance_chain, t_instance_chain = (
                query_system.protein_chains_asym_id[q_idx],
                target_protein_chains[t_idx],
            )
            q_chain, t_chain = (
                q_instance_chain.split(".")[1],
                t_instance_chain.split(".")[1],
            )
            try:
                aln = query_target_entry_alignments.loc[(q_chain, t_chain)]
            except KeyError:
                break
            pair = (q_instance_chain, t_instance_chain)
            pair_scores = self.get_protein_scores_pair(aln)
            for score in pair_scores:
                scores[f"{score}_weighted_sum"] += (
                    pair_scores[score]
                    * self.entries[query_system.pdb_id].chains[q_chain].length
                )
            alns[pair] = aln

            mappings[f"{protein_chain_mapper}_weighted_sum"].append(pair)
            # Zero out the entire row and column for the max element
            s_matrix[q_idx, :] = 0
            s_matrix[:, t_idx] = 0
        for score in scores:
            if score.endswith("_weighted_sum") and query_system_length > 0:
                scores[score] /= query_system_length
        return mappings, scores, alns, protein_chain_mapper

    def get_pocket_pli_scores(
        self,
        alns: dict[_ChainPairType, pd.DataFrame],
        query_system: System,
        target_system: System | None = None,
    ) -> tuple[
        _SimilarityScoreDictType,
        _SimilarityScoreDictType,
    ]:
        pocket_scores: _SimilarityScoreDictType = defaultdict(float)
        pli_scores: _SimilarityScoreDictType = defaultdict(float)
        pocket_length = query_system.num_pocket_residues
        pli_length = query_system.num_interactions
        pli_unique_length = query_system.num_unique_interactions
        for q_instance_chain, t_instance_chain in alns:
            aln = alns[(q_instance_chain, t_instance_chain)]
            q_pocket = query_system.pocket_residues.get(q_instance_chain, {})
            q_interactions = query_system.interactions_counter.get(q_instance_chain, {})
            if target_system is not None:
                t_pocket = target_system.pocket_residues.get(t_instance_chain, {})
                t_interactions = target_system.interactions_counter.get(
                    t_instance_chain, {}
                )
            for source, aln_source in aln.iterrows():
                for i in aln_source["qrnum"]:
                    q_a, t_a, lddt, q_n, t_n = (
                        aln_source["qaln"][i],
                        aln_source["taln"][i],
                        aln_source["lddtfull"].get(i, 0),
                        aln_source["qrnum"][i],
                        aln_source["trnum"].get(i, None),
                    )
                    assert q_a != "-" and t_a != "-"
                    if q_n in q_pocket:
                        pocket_scores[f"pocket_lddt_{source}"] += lddt
                        if q_a == t_a:
                            pocket_scores[f"pocket_fident_{source}"] += 1
                        if (
                            target_system is not None
                            and t_n is not None
                            and t_n in t_pocket
                        ):
                            pocket_scores[f"pocket_qcov_{source}"] += 1
                            pocket_scores[f"pocket_lddt_qcov_{source}"] += lddt
                            if q_a == t_a:
                                pocket_scores[f"pocket_fident_qcov_{source}"] += 1
                            if q_n in q_interactions and t_n in t_interactions:
                                pli_scores[f"pli_qcov_{source}"] += sum(
                                    (q_interactions[q_n] & t_interactions[t_n]).values()
                                )
                                pli_scores[f"pli_unique_qcov_{source}"] += len(
                                    (
                                        set(q_interactions[q_n].values())
                                        & set(t_interactions[t_n].values())
                                    )
                                )
        for score in pocket_scores:
            pocket_scores[score] /= pocket_length
        for score in pli_scores:
            if "unique" in score:
                pli_scores[score] /= pli_unique_length
            else:
                pli_scores[score] /= pli_length
        return pocket_scores, pli_scores

    def get_scores(
        self, search_db: str, query_system: System, query_entry_alignments: pd.DataFrame
    ) -> abc.Generator[dict[str, str | float], None, None]:
        if search_db == "holo":
            return self.get_scores_holo(query_system, query_entry_alignments)
        elif search_db == "apo" or search_db == "pred":
            return self.get_scores_apo_pred(query_system, query_entry_alignments)
        else:
            raise ValueError(f"Invalid search_db: {search_db}")

    def get_scores_holo(
        self, query_system: System, query_entry_alignments: pd.DataFrame
    ) -> abc.Generator[dict[str, str | float], None, None]:
        query_system_length = sum(
            self.entries[query_system.pdb_id]
            .chains[q_instance_chain.split(".")[1]]
            .length
            for q_instance_chain in query_system.protein_chains_asym_id
        )
        query_instance_chains = [
            chain.split(".")[1] for chain in query_system.protein_chains_asym_id
        ]
        for target_entry in query_entry_alignments.index.get_level_values(
            "target_entry"
        ).unique():
            if target_entry not in self.entries:
                # No data for this target entry
                continue
            query_target_entry_alignments = query_entry_alignments.loc[target_entry]
            query_chain_mapped_values = set(
                query_target_entry_alignments.index.get_level_values(
                    "query_chain_mapped"
                )
            )
            if all(
                chain not in query_chain_mapped_values
                for chain in query_instance_chains
            ):
                # No alignments for this query system
                continue
            all_target_chains: set[str] = set()
            for q_chain in query_instance_chains:
                if q_chain in query_chain_mapped_values:
                    all_target_chains.update(
                        query_target_entry_alignments.loc[
                            q_chain
                        ].index.get_level_values("target_chain_mapped")
                    )

            for target_system_id in self.entries[target_entry].systems:
                target_system = self.entries[target_entry].systems[target_system_id]
                if target_system.system_type != "holo":
                    continue
                if target_system_id == query_system.id or all(
                    target_instance_chain.split(".")[1] not in all_target_chains
                    for target_instance_chain in target_system.protein_chains_asym_id
                ):
                    # Same as query system or No alignments for this target system
                    continue
                q_t_scores: dict[str, float] = {}

                # Protein score calculation
                (
                    q_t_mappings,
                    protein_scores,
                    alns,
                    protein_chain_mapper,
                ) = self.get_protein_scores(
                    query_target_entry_alignments,
                    query_system,
                    target_system.protein_chains_asym_id,
                    query_system_length,
                )
                if not len(protein_scores):
                    continue
                q_t_scores.update(protein_scores)

                # Pocket and PLI score calculation
                (
                    pocket_scores,
                    pli_scores,
                ) = self.get_pocket_pli_scores(alns, query_system, target_system)
                q_t_scores.update(pocket_scores)
                q_t_scores.update(pli_scores)
                q_t_scores_combined = combine_scores(
                    q_t_scores, q_t_mappings, protein_chain_mapper
                )
                q_t_scores_combined["target_system"] = target_system.id
                q_t_scores_combined["query_system"] = query_system.id
                yield q_t_scores_combined

    def _get_suffixes(self, s: str) -> list[str]:
        suffixes = ["_weighted_sum", "_weighted_max", "_max"]
        if not s.startswith("protein"):
            suffixes = [""]
        return suffixes

    def get_column_mapr(self) -> dict[str, list[str]]:
        mapr = {"info": [x for x in INFO_COLUMNS]}
        for s in SCORE_NAMES:
            suffixes = self._get_suffixes(s)
            for suffix in suffixes:
                mapr[f"{s}{suffix}"] = [
                    f"{s}{suffix}",
                    f"{s}{suffix}_source",
                    f"{s}{suffix}_mapping",
                ]
        return mapr

    def aggregate_scores(
        self,
        pdb_id: str,
        search_db: str = "holo",
    ) -> Optional[pd.DataFrame]:
        alignments = self.load_alignments(
            search_db=search_db,
            source_to_aln_file={
                f"{search_db}_{aln_type}": self.db_dir
                / f"{search_db}_{aln_type}"
                / "mapped_aln"
                / f"{pdb_id}.parquet"
                for aln_type in ["foldseek", "mmseqs"]
            },
        )
        if alignments.empty:
            return None
        column_mapr = self.get_column_mapr()
        pdb_vals = []
        for system in self.entries[pdb_id].systems.values():
            if system.system_type != "holo":
                continue
            for score_dict in self.get_scores(
                search_db, system, alignments.loc[pdb_id]
            ):
                info = {}
                grouped: dict[str, dict[str, str | float]] = {}
                # iterate over groups of columns instead of raw list of columns
                for metric, columns in column_mapr.items():
                    for column in columns:
                        score = score_dict.get(column)
                        if score is None:
                            continue
                        if metric == "info":
                            info[column] = score
                        else:
                            grouped.setdefault(metric, {})
                            abbr_column = column.replace(metric, "", 1).lstrip("_")
                            if not abbr_column:
                                abbr_column = "similarity"
                                if np.isnan(score):
                                    continue
                            grouped[metric][abbr_column] = score
                for metric, values in grouped.items():
                    if len(values) and values.get("similarity") is not None:
                        pdb_vals.append({**info, **values, **{"metric": metric}})
        if not len(pdb_vals):
            return None
        df = pd.DataFrame(pdb_vals)
        all_thresholds = [
            (m, self.minimum_thresholds.get(m, self.minimum_threshold))
            for m in df["metric"].unique()
        ]
        query = " or ".join(
            [f"(metric=='{m}' and similarity>={t})" for (m, t) in all_thresholds]
        )
        df = df.query(query).copy()
        for col in [
            "protein_mapper",
            "source",
            "metric",
        ]:
            df[col] = df[col].astype("category")
        df["similarity"] = (df["similarity"] * 100).apply(round).astype(np.int8)
        columns = [col for (col, _) in SORT_ORDER]
        ascending = [direction == "ascending" for (_, direction) in SORT_ORDER]
        return df.sort_values(by=columns, ascending=ascending)

    def get_scores_apo_pred(
        self, query_system: System, query_entry_alignments: pd.DataFrame
    ) -> abc.Generator[dict[str, str | float], None, None]:
        q_chain = query_system.protein_chains_asym_id[0].split(".")[1]
        query_system_length = sum(
            self.entries[query_system.pdb_id]
            .chains[q_instance_chain.split(".")[1]]
            .length
            for q_instance_chain in query_system.protein_chains_asym_id
        )
        for target_entry in query_entry_alignments.index.get_level_values(
            "target_entry"
        ).unique():
            query_target_entry_alignments = query_entry_alignments.loc[target_entry]
            try:
                q_chain_alignments = query_target_entry_alignments.loc[q_chain]
            except KeyError:
                q_chain_alignments = None
            if q_chain_alignments is not None:
                for t_chain in q_chain_alignments.index.get_level_values(
                    "target_chain_mapped"
                ).unique():
                    if (
                        target_entry in self.entries
                        and self.entries[target_entry].chains[t_chain].holo
                    ):
                        continue
                    q_t_scores: dict[str, float] = {}
                    # Protein score calculation
                    (
                        q_t_mappings,
                        protein_scores,
                        alns,
                        protein_chain_mapper,
                    ) = self.get_protein_scores(
                        query_target_entry_alignments,
                        query_system,
                        [f"0.{t_chain}"],
                        query_system_length,
                    )

                    if len(alns) == 0 or not len(protein_scores):
                        continue
                    q_t_scores.update(protein_scores)
                    # Pocket score calculation
                    q_t_scores.update(self.get_pocket_pli_scores(alns, query_system)[0])
                    q_t_scores_combined = combine_scores(
                        q_t_scores, q_t_mappings, protein_chain_mapper
                    )
                    q_t_scores_combined["target_system"] = f"{target_entry}_{t_chain}"
                    q_t_scores_combined["query_system"] = query_system.id
                    yield q_t_scores_combined
