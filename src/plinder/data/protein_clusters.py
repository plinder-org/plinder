# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
"""Native protein-chain clustering for release builds."""

from __future__ import annotations

import hashlib
import json
import subprocess
from concurrent.futures import ThreadPoolExecutor
from itertools import islice
from pathlib import Path
from tempfile import TemporaryDirectory

import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

from plinder.core.release import RELEASE_PATHS
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


def _clusters_are_current(output: Path, metadata: bytes) -> bool:
    if not output.is_file():
        return False
    return (pq.read_schema(output).metadata or {}).get(CLUSTER_METADATA_KEY) == metadata


def _read_assignments(path: Path, expected_ids: set[str]) -> pd.DataFrame:
    """Require one native-tool assignment per submitted input."""
    assignments = pd.read_csv(
        path,
        sep="\t",
        header=None,
        names=["representative", "member"],
        dtype=str,
        keep_default_na=False,
    )
    observed = set(assignments["member"])
    if len(assignments) != len(expected_ids) or observed != expected_ids:
        raise ValueError(
            "protein cluster membership differs from input: "
            f"rows={len(assignments)}, expected={len(expected_ids)}, "
            f"missing={sorted(expected_ids - observed)[:10]}, "
            f"extra={sorted(observed - expected_ids)[:10]}"
        )
    return assignments


def _write_cluster_table(
    chains: pd.DataFrame,
    assignments: pd.DataFrame,
    output: Path,
    metadata: bytes,
    *,
    missing_status: str | dict[str, str],
) -> None:
    representative_keys = chains[["member", "entry_pdb_id", "chain_asym_id"]].rename(
        columns={
            "member": "representative",
            "entry_pdb_id": "representative_entry_pdb_id",
            "chain_asym_id": "representative_chain_asym_id",
        }
    )
    result = chains.merge(assignments, on="member", how="left", validate="one_to_one")
    result = result.merge(
        representative_keys, on="representative", how="left", validate="many_to_one"
    )
    result["is_representative"] = (
        result["member"].eq(result["representative"]).astype("boolean")
    )
    missing = result["representative"].isna()
    result.loc[missing, "is_representative"] = pd.NA
    result["status"] = "clustered"
    result.loc[missing, "status"] = (
        result.loc[missing, "member"].map(missing_status)
        if isinstance(missing_status, dict)
        else missing_status
    )
    table = pa.Table.from_pandas(
        result[SEQUENCE_CLUSTER_SCHEMA.names],
        schema=SEQUENCE_CLUSTER_SCHEMA,
        preserve_index=False,
    ).replace_schema_metadata({CLUSTER_METADATA_KEY: metadata})
    output.parent.mkdir(parents=True, exist_ok=True)
    # Stage beside the destination so replacement is atomic across filesystems.
    with TemporaryDirectory(prefix=f".{output.stem}-", dir=output.parent) as staging:
        staged = Path(staging) / output.name
        pq.write_table(table, staged, compression="zstd")
        staged.replace(output)


def _prepare_structure_entry(
    cif_file: Path, chains: pd.DataFrame, input_dir: Path
) -> dict[str, str]:
    """Write usable label-asym chains and return each chain's preparation status."""
    import numpy as np

    from plinder.data.annotations.cif_utils import (
        get_structure_with_altloc,
        read_mmcif_file,
    )
    from plinder.data.annotations.save_utils import save_cif_file

    source = read_mmcif_file(cif_file)
    atoms = get_structure_with_altloc(source)
    statuses = {}
    for row in chains.itertuples(index=False):
        selected = atoms[atoms.chain_id == row.chain_asym_id].copy()
        ca = selected[(selected.atom_name == "CA") & (selected.element == "C")]
        # A 3Di description needs a local backbone neighbourhood.
        if len(ca) < 4:
            statuses[row.member] = "insufficient_coordinates"
            continue
        if not np.isfinite(selected.coord).all():
            raise ValueError(
                f"non-finite coordinates in {cif_file}, chain {row.chain_asym_id}"
            )
        # Foldseek rejects chains whose resolved residues are all unknown.
        if np.all(ca.res_name == "UNK"):
            statuses[row.member] = "unknown_residues"
            continue
        selected.chain_id[:] = "A"
        # With --chain-name-mode 1, Foldseek appends the internal chain ID.
        filename = row.member.removesuffix("_A")
        save_cif_file(
            selected,
            filename,
            input_dir / f"{filename}.cif",
            source_block=source.block,
            source_asym_ids={"A": row.chain_asym_id},
        )
        statuses[row.member] = "clustered"
    return statuses


def make_protein_structure_clusters(
    *,
    data_dir: Path,
    cif_root: Path,
    scratch_dir: Path,
    threads: int = 4,
    lddt: float = 0.7,
    coverage: float = 0.8,
    force_update: bool = False,
) -> Path:
    """Cluster release protein chains with Foldseek easy-cluster.

    Coverage applies to both resolved chains. The first model and deposited-first
    alternate conformers match entry ingest. Chains with fewer than four resolved
    C-alpha atoms retain a null assignment and ``insufficient_coordinates`` status.
    Chains whose resolved residues are all UNK have ``unknown_residues`` status.
    Missing source files and parsing failures stop the stage. Inputs and native
    clustering intermediates are temporary; only the chain table is retained.
    """
    from plinder.data.pipeline.ingest import resolve_entry_paths

    if threads < 1:
        raise ValueError("protein clustering threads must be positive")
    if not 0 <= lddt <= 1 or not 0 <= coverage <= 1:
        raise ValueError("protein clustering lDDT and coverage must be between 0 and 1")
    chains = (
        pd.read_parquet(
            data_dir / "index/entry_chains.parquet",
            columns=["entry_pdb_id", "chain_asym_id"],
            filters=[("chain_receptor_type", "==", "protein")],
        )
        .sort_values(["entry_pdb_id", "chain_asym_id"])
        .reset_index(drop=True)
    )
    keys = ["entry_pdb_id", "chain_asym_id"]
    chains["member"] = pd.Series(
        [f"chain_{i}_A" for i in range(len(chains))], dtype="string"
    )
    digest = hashlib.sha256()
    digest.update(chains[keys].to_json(orient="values").encode())
    sources = {}
    for pdb_id in chains["entry_pdb_id"].unique():
        path, _ = resolve_entry_paths(
            pdb_id, cif_root=cif_root, validation_root=cif_root
        )
        stat = path.stat()
        sources[pdb_id] = path
        digest.update(
            json.dumps([str(path.resolve()), stat.st_size, stat.st_mtime_ns]).encode()
        )
        digest.update(b"\n")
    parameters = {
        "backend": "foldseek",
        "version": subprocess.check_output(["foldseek", "version"], text=True).strip(),
        "lddt": lddt,
        "coverage": coverage,
        "coverage_mode": 0,
        "cluster_mode": 0,
        "single_step_clustering": True,
        "alignment_type": 2,
        "chain_name_mode": 1,
        "model": 1,
        "altloc": "first",
        "minimum_ca_atoms": 4,
        "exclude_unknown_only_chains": True,
        "source_signature": digest.hexdigest(),
    }
    metadata = json.dumps(parameters, sort_keys=True).encode()
    output = data_dir / RELEASE_PATHS["protein_structure_clusters"]
    if not force_update and _clusters_are_current(output, metadata):
        return output
    scratch_dir.mkdir(parents=True, exist_ok=True)
    with TemporaryDirectory(
        prefix="plinder-structure-clusters-", dir=scratch_dir
    ) as temporary:
        work = Path(temporary)
        input_dir = work / "chains"
        input_dir.mkdir()
        statuses: dict[str, str] = {}
        groups = iter(chains.groupby("entry_pdb_id", sort=False))
        with ThreadPoolExecutor(max_workers=threads) as executor:
            # Bound queued entries and memory use, including on Python < 3.14.
            while batch := list(islice(groups, threads)):
                futures = [
                    executor.submit(
                        _prepare_structure_entry, sources[pdb_id], group, input_dir
                    )
                    for pdb_id, group in batch
                ]
                for future in futures:
                    statuses.update(future.result())
        prepared = {
            member for member, status in statuses.items() if status == "clustered"
        }
        if prepared:
            run(
                [
                    "foldseek",
                    "easy-cluster",
                    str(input_dir),
                    str(work / "clusters"),
                    str(work / "tmp"),
                    "--lddt-threshold",
                    str(lddt),
                    "-c",
                    str(coverage),
                    "--cov-mode",
                    "0",
                    "--cluster-mode",
                    "0",
                    "--single-step-clustering",
                    "1",
                    "--alignment-type",
                    "2",
                    "--chain-name-mode",
                    "1",
                    "--threads",
                    str(threads),
                ]
            )
            assignments = _read_assignments(work / "clusters_cluster.tsv", prepared)
        else:
            assignments = pd.DataFrame(columns=["representative", "member"])
        _write_cluster_table(
            chains,
            assignments,
            output,
            metadata,
            missing_status=statuses,
        )
    return output


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
    output = data_dir / RELEASE_PATHS["protein_sequence_clusters"]
    if not force_update and _clusters_are_current(output, metadata):
        return output
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
        _write_cluster_table(
            chains, assignments, output, metadata, missing_status="missing_sequence"
        )
    return output
