# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
"""Native protein-chain clustering for release builds."""

from __future__ import annotations

import hashlib
import json
import subprocess
from pathlib import Path
from tempfile import TemporaryDirectory

import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

from plinder.data.databases import run

SEQUENCE_CLUSTER_SCHEMA = pa.schema(
    [
        ("entry_pdb_id", pa.string()),
        ("chain_asym_id", pa.string()),
        ("representative_entry_pdb_id", pa.string()),
        ("representative_chain_asym_id", pa.string()),
        ("is_representative", pa.bool_()),
        ("status", pa.string()),
    ]
)
CLUSTER_METADATA_KEY = b"plinder_protein_clustering"


def _read_assignments(path: Path, expected_ids: set[str]) -> pd.DataFrame:
    """Check that every input has one assignment to a retained representative."""
    assignments = pd.read_csv(
        path,
        sep="\t",
        header=None,
        names=["representative", "member"],
        dtype=str,
        keep_default_na=False,
    )
    if assignments.empty or assignments["member"].duplicated().any():
        raise ValueError("protein clusters must assign each input exactly once")
    observed = set(assignments["member"])
    if observed != expected_ids:
        raise ValueError(
            "protein cluster membership differs from input: "
            f"missing={sorted(expected_ids - observed)[:10]}, "
            f"extra={sorted(observed - expected_ids)[:10]}"
        )
    representatives = set(assignments["representative"])
    self_members = set(
        assignments.loc[
            assignments["representative"] == assignments["member"], "member"
        ]
    )
    if representatives != self_members:
        raise ValueError(
            "protein cluster representatives must belong to their own clusters"
        )
    return assignments


def make_protein_sequence_clusters(
    *,
    data_dir: Path,
    scratch_dir: Path,
    threads: int = 4,
    identity: float = 0.4,
    coverage: float = 0.8,
    force_update: bool = False,
) -> Path:
    """Cluster all release protein-chain sequences using MMseqs easy-cluster.

    Each representative/member match must meet the identity threshold and cover
    the requested fraction of both chains. Greedy set cover runs in single-step
    mode. Apo and interface chains are included. Missing sequences have null
    assignments and ``status='missing_sequence'`` in the output.

    Intermediate MMseqs files live under ``scratch_dir`` and are removed after
    the table is written. Parameters, input sequence digest, and MMseqs version
    are stored in Parquet metadata to check whether a previous result is reusable.
    """
    if threads < 1:
        raise ValueError("protein clustering threads must be positive")
    if not 0 <= identity <= 1 or not 0 <= coverage <= 1:
        raise ValueError(
            "protein clustering identity and coverage must be between 0 and 1"
        )
    chains = (
        pd.read_parquet(
            data_dir / "index" / "entry_chains.parquet",
            columns=["entry_pdb_id", "chain_asym_id", "chain_sequence"],
            filters=[("chain_receptor_type", "==", "protein")],
        )
        .sort_values(["entry_pdb_id", "chain_asym_id"])
        .reset_index(drop=True)
    )
    keys = ["entry_pdb_id", "chain_asym_id"]
    if chains[keys].isna().any().any() or chains.duplicated(keys).any():
        raise ValueError("protein chains require unique, non-null entry/asym IDs")
    # Generated FASTA IDs avoid any restrictions in native tools on source IDs.
    # They are mapped back to the untouched entry/asym IDs before writing.
    chains["member"] = pd.Series(
        [f"chain_{i}" for i in range(len(chains))], dtype="string"
    )
    valid = chains["chain_sequence"].notna() & chains["chain_sequence"].ne("")
    sequences = chains.loc[valid]
    if not sequences["chain_sequence"].str.fullmatch("[A-Za-z]+").all():
        raise ValueError(
            "protein chain sequences must contain one-letter amino-acid codes"
        )
    digest = hashlib.sha256()
    for row in chains[keys + ["chain_sequence"]].itertuples(index=False, name=None):
        digest.update(
            json.dumps([None if pd.isna(value) else value for value in row]).encode()
        )
        digest.update(b"\n")
    parameters = {
        "backend": "mmseqs",
        "sequence_type": "protein",
        "version": subprocess.check_output(["mmseqs", "version"], text=True).strip(),
        "identity": identity,
        "coverage": coverage,
        "coverage_mode": 0,
        "cluster_mode": 0,
        "single_step_clustering": True,
        "sequence_digest": digest.hexdigest(),
    }
    metadata = json.dumps(parameters, sort_keys=True).encode()
    output = data_dir / "protein_clusters" / "sequence.parquet"
    if not force_update and output.is_file():
        try:
            previous = pq.read_table(output)
            if (
                previous.schema.equals(SEQUENCE_CLUSTER_SCHEMA)
                and previous.num_rows == len(chains)
                and (previous.schema.metadata or {}).get(CLUSTER_METADATA_KEY)
                == metadata
            ):
                return output
        except (OSError, ValueError, pa.ArrowException):
            pass
    scratch_dir.mkdir(parents=True, exist_ok=True)
    with TemporaryDirectory(
        prefix="plinder-sequence-clusters-", dir=scratch_dir
    ) as temporary:
        work = Path(temporary)
        if sequences.empty:
            assignments = pd.DataFrame(columns=["representative", "member"])
        else:
            fasta = work / "chains.fasta"
            with fasta.open("w") as handle:
                for row in sequences.itertuples(index=False):
                    handle.write(f">{row.member}\n{row.chain_sequence}\n")
            run(
                [
                    "mmseqs",
                    "easy-cluster",
                    str(fasta),
                    str(work / "clusters"),
                    str(work / "tmp"),
                    "--dbtype",
                    "1",
                    "--min-seq-id",
                    str(identity),
                    "-c",
                    str(coverage),
                    "--cov-mode",
                    "0",
                    "--cluster-mode",
                    "0",
                    "--single-step-clustering",
                    "1",
                    "--threads",
                    str(threads),
                ]
            )
            assignments = _read_assignments(
                work / "clusters_cluster.tsv", set(sequences["member"])
            )
        representative_keys = chains[["member", *keys]].rename(
            columns={
                "member": "representative",
                "entry_pdb_id": "representative_entry_pdb_id",
                "chain_asym_id": "representative_chain_asym_id",
            }
        )
        result = chains.merge(
            assignments, on="member", how="left", validate="one_to_one"
        )
        result = result.merge(
            representative_keys, on="representative", how="left", validate="many_to_one"
        )
        result["is_representative"] = (
            result["member"].eq(result["representative"]).astype("boolean")
        )
        missing = result["representative"].isna()
        result.loc[missing, "is_representative"] = pd.NA
        result["status"] = "clustered"
        result.loc[missing, "status"] = "missing_sequence"
        table = pa.Table.from_pandas(
            result[SEQUENCE_CLUSTER_SCHEMA.names],
            schema=SEQUENCE_CLUSTER_SCHEMA,
            preserve_index=False,
        ).replace_schema_metadata({CLUSTER_METADATA_KEY: metadata})
        output.parent.mkdir(parents=True, exist_ok=True)
        # Stage beside the destination so replacement is atomic across filesystems.
        with TemporaryDirectory(prefix=".sequence-", dir=output.parent) as staging:
            staged = Path(staging) / output.name
            pq.write_table(table, staged, compression="zstd")
            staged.replace(output)
    return output
