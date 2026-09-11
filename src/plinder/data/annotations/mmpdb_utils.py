# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
"""Build a compact matched-molecular-pair table for unique ligand SMILES."""

from __future__ import annotations

import gzip
import hashlib
import json
import shutil
import subprocess
from concurrent.futures import ThreadPoolExecutor
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Sequence

import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
from rdkit import Chem
from rdkit.rdBase import BlockLogs

from plinder.core.utils import schemas
from plinder.core.utils.log import setup_logger

LOG = setup_logger(__name__)


def _mmpdb_version() -> str:
    try:
        return version("mmpdb")
    except PackageNotFoundError:
        return "unknown"


def _empty_pair_table() -> pa.Table:
    return pa.Table.from_arrays(
        [pa.array([], type=field.type) for field in schemas.LIGAND_MMP_PAIR_SCHEMA],
        schema=schemas.LIGAND_MMP_PAIR_SCHEMA,
    )


def _ligand_table(fingerprint_path: Path) -> pd.DataFrame:
    if not fingerprint_path.is_file():
        raise FileNotFoundError(
            "matched molecular pairs require the unique ligand fingerprint table: "
            f"{fingerprint_path}"
        )
    ligands = pd.read_parquet(
        fingerprint_path,
        columns=["ligand_smiles_id", "ligand_rdkit_canonical_smiles"],
    )
    if ligands["ligand_smiles_id"].isna().any():
        raise ValueError("unique ligand SMILES contain missing ligand_smiles_id values")
    if ligands["ligand_rdkit_canonical_smiles"].isna().any():
        raise ValueError("unique ligand SMILES contain missing structures")
    if ligands["ligand_smiles_id"].duplicated().any():
        raise ValueError(
            "unique ligand SMILES contain duplicate ligand_smiles_id values"
        )
    if ligands["ligand_rdkit_canonical_smiles"].duplicated().any():
        raise ValueError("unique ligand SMILES contain duplicate structures")

    ligands = ligands.sort_values("ligand_smiles_id").reset_index(drop=True)
    ligands["ligand_smiles_id"] = ligands["ligand_smiles_id"].astype("int32")
    ligands["ligand_rdkit_canonical_smiles"] = ligands[
        "ligand_rdkit_canonical_smiles"
    ].astype(str)
    heavy_atom_counts: list[int] = []
    contains_dative_bond: list[bool] = []
    with BlockLogs():
        for ligand_id, smiles in ligands[
            ["ligand_smiles_id", "ligand_rdkit_canonical_smiles"]
        ].itertuples(index=False, name=None):
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                raise ValueError(
                    f"ligand_smiles_id {ligand_id} has invalid canonical SMILES"
                )
            heavy_atom_counts.append(mol.GetNumHeavyAtoms())
            contains_dative_bond.append(
                any(
                    str(bond.GetBondType()).startswith("DATIVE")
                    for bond in mol.GetBonds()
                )
            )
    ligands["num_heavy_atoms"] = pd.Series(heavy_atom_counts, dtype="int16")
    if any(contains_dative_bond):
        LOG.info(
            "excluding %d ligand SMILES with metal-dative bonds from MMP generation",
            sum(contains_dative_bond),
        )
        ligands = ligands.loc[
            ~pd.Series(contains_dative_bond, index=ligands.index)
        ].reset_index(drop=True)
    return ligands


def _ligand_signature(ligands: pd.DataFrame) -> str:
    digest = hashlib.sha256()
    for ligand_id, smiles in ligands[
        ["ligand_smiles_id", "ligand_rdkit_canonical_smiles"]
    ].itertuples(index=False, name=None):
        digest.update(str(int(ligand_id)).encode())
        digest.update(b"\0")
        digest.update(smiles.encode())
        digest.update(b"\n")
    return digest.hexdigest()


def _run_command(arguments: Sequence[str], *, cwd: Path) -> None:
    LOG.info("running %s", " ".join(arguments))
    try:
        completed = subprocess.run(
            list(arguments),
            cwd=cwd,
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
        )
    except subprocess.CalledProcessError as exc:
        output = (exc.stdout or "").strip()
        if len(output) > 8_000:
            output = output[-8_000:]
        raise RuntimeError(
            f"mmpdb command failed ({' '.join(arguments)}):\n{output}"
        ) from exc
    if completed.stdout:
        LOG.info("mmpdb: %s", completed.stdout.rstrip())


def _run_commands(
    commands: Sequence[Sequence[str]], *, cwd: Path, workers: int
) -> None:
    if not commands:
        return
    with ThreadPoolExecutor(max_workers=min(workers, len(commands))) as executor:
        futures = [
            executor.submit(_run_command, command, cwd=cwd) for command in commands
        ]
        for future in futures:
            future.result()


def _gzip_has_rows(path: Path) -> bool:
    with gzip.open(path, "rb") as handle:
        return bool(handle.read(1))


def _generate_pair_files(
    *, ligands: pd.DataFrame, work_dir: Path, threads: int, executable: str
) -> list[Path]:
    smiles_path = work_dir / "ligands.smi"
    ligands[["ligand_rdkit_canonical_smiles", "ligand_smiles_id"]].to_csv(
        smiles_path, sep="\t", header=False, index=False
    )
    if ligands.empty:
        return []

    shard_count = min(threads, len(ligands))
    _run_command(
        [executable, "smi_split", "-n", str(shard_count), smiles_path.name],
        cwd=work_dir,
    )
    smiles_shards = sorted(work_dir.glob("ligands.[0-9][0-9][0-9][0-9].smi"))
    if not smiles_shards:
        raise RuntimeError("mmpdb did not create any SMILES shards")
    _run_commands(
        [[executable, "fragment", "-j", "1", shard.name] for shard in smiles_shards],
        cwd=work_dir,
        workers=threads,
    )
    fragment_files = sorted(work_dir.glob("ligands.[0-9][0-9][0-9][0-9].fragdb"))
    if len(fragment_files) != len(smiles_shards):
        raise RuntimeError(
            "mmpdb fragmentation did not produce one database per SMILES shard"
        )

    _run_command(
        [
            executable,
            "fragdb_partition",
            "-n",
            str(shard_count),
            "--template",
            "partition.{i:04}.fragdb",
            *[path.name for path in fragment_files],
        ],
        cwd=work_dir,
    )
    partition_files = sorted(work_dir.glob("partition.[0-9][0-9][0-9][0-9].fragdb"))
    if not partition_files:
        return []
    _run_commands(
        [
            [
                executable,
                "index",
                "--out",
                "csv.gz",
                "-o",
                partition.with_suffix(".csv.gz").name,
                partition.name,
            ]
            for partition in partition_files
        ],
        cwd=work_dir,
        workers=threads,
    )
    pair_files = sorted(work_dir.glob("partition.[0-9][0-9][0-9][0-9].csv.gz"))
    if len(pair_files) != len(partition_files):
        raise RuntimeError("mmpdb indexing did not produce one pair file per partition")
    return [path for path in pair_files if _gzip_has_rows(path)]


def _core_table(shared_cores: Sequence[str]) -> pd.DataFrame:
    rows: list[tuple[str, int, int]] = []
    with BlockLogs():
        for shared_core in shared_cores:
            mol = Chem.MolFromSmiles(shared_core)
            if mol is None:
                mol = Chem.MolFromSmarts(shared_core)
            if mol is None:
                raise ValueError(
                    f"mmpdb emitted an unreadable shared core: {shared_core}"
                )
            rows.append((shared_core, shared_core.count("*"), mol.GetNumHeavyAtoms()))
    return pd.DataFrame(
        rows,
        columns=[
            "shared_core_smiles",
            "num_cuts",
            "shared_core_num_heavy_atoms",
        ],
    ).astype({
        "num_cuts": "int8",
        "shared_core_num_heavy_atoms": "int16",
    })


def _sql_paths(paths: Sequence[Path]) -> str:
    return ", ".join(
        f"'{path.as_posix().replace(chr(39), chr(39) * 2)}'" for path in paths
    )


def _write_pair_parquet(
    *, pair_files: Sequence[Path], ligands: pd.DataFrame, output_path: Path
) -> None:
    if not pair_files:
        pq.write_table(_empty_pair_table(), output_path, compression="zstd")
        return

    import duckdb

    connection = duckdb.connect()
    raw_path = pair_files[0].parent / "raw_pairs.parquet"
    raw_path_sql = raw_path.as_posix().replace("'", "''")
    pair_paths_sql = _sql_paths(pair_files)
    connection.execute(
        f"""
        COPY (
            SELECT
                try_cast(ligand_smiles_id_1 AS INTEGER) AS ligand_smiles_id_1,
                try_cast(ligand_smiles_id_2 AS INTEGER) AS ligand_smiles_id_2,
                transformation,
                shared_core_smiles
            FROM read_csv(
                [{pair_paths_sql}],
                delim='\t',
                header=false,
                columns={{
                    'raw_smiles_1': 'VARCHAR',
                    'raw_smiles_2': 'VARCHAR',
                    'ligand_smiles_id_1': 'VARCHAR',
                    'ligand_smiles_id_2': 'VARCHAR',
                    'transformation': 'VARCHAR',
                    'shared_core_smiles': 'VARCHAR'
                }}
            )
        ) TO '{raw_path_sql}'
        (FORMAT PARQUET, COMPRESSION ZSTD, ROW_GROUP_SIZE 1_000_000)
        """
    )
    invalid_ids = connection.execute(
        f"""
        SELECT count(*)
        FROM read_parquet('{raw_path_sql}')
        WHERE ligand_smiles_id_1 IS NULL OR ligand_smiles_id_2 IS NULL
        """
    ).fetchone()[0]
    if invalid_ids:
        connection.close()
        raise ValueError(f"mmpdb emitted {invalid_ids} non-integer ligand IDs")

    shared_cores = (
        connection
        .execute(
            f"""
        SELECT DISTINCT shared_core_smiles
        FROM read_parquet('{raw_path_sql}')
        ORDER BY shared_core_smiles
        """
        )
        .fetch_df()["shared_core_smiles"]
        .tolist()
    )
    core_table = _core_table(shared_cores)
    ligand_lookup = ligands.rename(
        columns={
            "ligand_rdkit_canonical_smiles": "ligand_smiles",
            "num_heavy_atoms": "ligand_num_heavy_atoms",
        }
    )[["ligand_smiles_id", "ligand_smiles", "ligand_num_heavy_atoms"]]
    connection.register("ligand_lookup", ligand_lookup)
    connection.register("core_lookup", core_table)
    missing_ids = connection.execute(
        f"""
        SELECT count(*)
        FROM read_parquet('{raw_path_sql}') raw
        LEFT JOIN ligand_lookup ligand_1
          ON raw.ligand_smiles_id_1 = ligand_1.ligand_smiles_id
        LEFT JOIN ligand_lookup ligand_2
          ON raw.ligand_smiles_id_2 = ligand_2.ligand_smiles_id
        WHERE ligand_1.ligand_smiles_id IS NULL
           OR ligand_2.ligand_smiles_id IS NULL
        """
    ).fetchone()[0]
    if missing_ids:
        connection.close()
        raise ValueError(
            f"mmpdb emitted {missing_ids} pairs with unknown ligand_smiles_id values"
        )

    output_sql = output_path.as_posix().replace("'", "''")
    connection.execute(
        f"""
        COPY (
            SELECT DISTINCT
                raw.ligand_smiles_id_1::INTEGER AS ligand_smiles_id_1,
                raw.ligand_smiles_id_2::INTEGER AS ligand_smiles_id_2,
                ligand_1.ligand_smiles::VARCHAR AS ligand_smiles_1,
                ligand_2.ligand_smiles::VARCHAR AS ligand_smiles_2,
                raw.transformation::VARCHAR AS transformation,
                raw.shared_core_smiles::VARCHAR AS shared_core_smiles,
                core.num_cuts::TINYINT AS num_cuts,
                core.shared_core_num_heavy_atoms::SMALLINT
                    AS shared_core_num_heavy_atoms,
                ligand_1.ligand_num_heavy_atoms::SMALLINT
                    AS ligand_1_num_heavy_atoms,
                ligand_2.ligand_num_heavy_atoms::SMALLINT
                    AS ligand_2_num_heavy_atoms,
                cast(
                    core.shared_core_num_heavy_atoms::FLOAT
                    / nullif(ligand_1.ligand_num_heavy_atoms, 0)
                    AS FLOAT
                ) AS ligand_1_shared_core_fraction,
                cast(
                    core.shared_core_num_heavy_atoms::FLOAT
                    / nullif(ligand_2.ligand_num_heavy_atoms, 0)
                    AS FLOAT
                ) AS ligand_2_shared_core_fraction
            FROM read_parquet('{raw_path_sql}') raw
            INNER JOIN ligand_lookup ligand_1
              ON raw.ligand_smiles_id_1 = ligand_1.ligand_smiles_id
            INNER JOIN ligand_lookup ligand_2
              ON raw.ligand_smiles_id_2 = ligand_2.ligand_smiles_id
            INNER JOIN core_lookup core USING (shared_core_smiles)
            ORDER BY
                ligand_smiles_id_1,
                ligand_smiles_id_2,
                shared_core_smiles,
                transformation
        ) TO '{output_sql}'
        (FORMAT PARQUET, COMPRESSION ZSTD, ROW_GROUP_SIZE 1_000_000)
        """
    )
    connection.close()

    actual_schema = pq.read_schema(output_path)
    if not actual_schema.equals(schemas.LIGAND_MMP_PAIR_SCHEMA):
        raise ValueError(
            f"ligand MMP output schema differs from the release schema: {actual_schema}"
        )


def make_ligand_mmp_pairs(
    *,
    data_dir: Path,
    scratch_dir: Path,
    threads: int = 4,
    force_update: bool = False,
) -> Path:
    """Generate matched molecular pairs for every unique proper-ligand SMILES."""
    if threads < 1:
        raise ValueError("MMP generation threads must be positive")
    fingerprint_path = data_dir / "fingerprints" / "ligands_per_smiles.parquet"
    output_path = data_dir / "index" / "ligand_mmp_pairs.parquet"
    manifest_path = data_dir / "index" / "ligand_mmp_pairs.manifest.json"
    ligands = _ligand_table(fingerprint_path)
    manifest = {
        "mmpdb_version": _mmpdb_version(),
        "ligand_signature": _ligand_signature(ligands),
    }
    if not force_update and output_path.is_file() and manifest_path.is_file():
        try:
            cached_manifest = json.loads(manifest_path.read_text())
            schema_is_current = pq.read_schema(output_path).equals(
                schemas.LIGAND_MMP_PAIR_SCHEMA
            )
        except (OSError, ValueError, json.JSONDecodeError):
            schema_is_current = False
            cached_manifest = None
        if cached_manifest == manifest and schema_is_current:
            LOG.info("ligand MMP pairs already match the unique ligand SMILES")
            return output_path

    executable = shutil.which("mmpdb")
    if executable is None:
        raise RuntimeError(
            "mmpdb is required to generate ligand MMP pairs; install the PLINDER "
            "data dependencies"
        )
    output_path.parent.mkdir(exist_ok=True, parents=True)
    scratch_dir.mkdir(exist_ok=True, parents=True)
    temporary_output = output_path.with_suffix(".tmp.parquet")
    temporary_output.unlink(missing_ok=True)
    try:
        with TemporaryDirectory(prefix="plinder-mmp-", dir=scratch_dir) as work:
            pair_files = _generate_pair_files(
                ligands=ligands,
                work_dir=Path(work),
                threads=threads,
                executable=executable,
            )
            _write_pair_parquet(
                pair_files=pair_files,
                ligands=ligands,
                output_path=temporary_output,
            )
        temporary_output.replace(output_path)
    except BaseException:
        temporary_output.unlink(missing_ok=True)
        raise

    temporary_manifest = manifest_path.with_suffix(".tmp.json")
    temporary_manifest.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    temporary_manifest.replace(manifest_path)
    LOG.info("wrote ligand MMP pairs to %s", output_path)
    return output_path
