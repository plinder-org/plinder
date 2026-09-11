# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
"""
Wrap all network requests in a retry decorator
and use a convention of looking for a file in a
pre-determined location before fetching it from
the network.
"""

import gzip
import re
from collections.abc import Iterable
from concurrent.futures import ThreadPoolExecutor
from functools import wraps
from pathlib import Path
from tempfile import NamedTemporaryFile
from time import sleep
from typing import Any, Callable, Optional, TypeVar

import pandas as pd
import requests
from biotite.database.rcsb import fetch
from biotite.file import DeserializationError, InvalidFileError
from biotite.structure.io.pdbx import CIFFile, get_structure, set_structure

from plinder.core.release import PlinderRelease
from plinder.core.utils.config import get_config
from plinder.core.utils.cpl import is_offline
from plinder.core.utils.log import setup_logger

LOG = setup_logger(__name__)
T = TypeVar("T")
WWPDB_VERSIONED_MMCIF_URL = (
    "https://files-versioned.wwpdb.org/pdb_versioned/data/entries"
)
PDBRevision = tuple[int, int]
SOURCE_MMCIF_DOWNLOAD_ATTEMPTS = 3
TRANSIENT_HTTP_STATUSES = frozenset({429, *range(500, 600)})


def _normalize_pdb_id(pdb_or_system_id: str) -> str:
    """Return a safe lowercase PDB identifier from a PDB or system ID."""
    pdb_id = str(pdb_or_system_id).split("__", maxsplit=1)[0].strip().lower()
    if not re.fullmatch(r"(?:[0-9][a-z0-9]{3}|pdb_[a-z0-9]{8})", pdb_id):
        raise ValueError(f"invalid PDB ID: {pdb_or_system_id!r}")
    # The extended representation of a legacy four-character accession is
    # equivalent to its short form.  Canonicalizing it avoids duplicate cache
    # entries and manifest lookups.
    if pdb_id.startswith("pdb_0000"):
        return pdb_id[-4:]
    return pdb_id


def _entry_source_manifest_path(manifest_path: Path | str | None) -> Path:
    if manifest_path is not None:
        return Path(manifest_path)
    return PlinderRelease().fetch("entry_sources")


def _get_pdb_revisions(
    pdb_ids: Iterable[str],
    *,
    manifest_path: Path | str | None = None,
) -> dict[str, PDBRevision]:
    """Read exact source revisions for a subset from the release manifest."""
    normalized = list(dict.fromkeys(_normalize_pdb_id(value) for value in pdb_ids))
    if not normalized:
        return {}
    path = _entry_source_manifest_path(manifest_path)
    if not path.is_file():
        raise FileNotFoundError(
            f"source mmCIF revision manifest is not available at {path}; "
            "download the release index before entering offline mode"
        )
    columns = [
        "entry_pdb_id",
        "source_mmcif_major_revision",
        "source_mmcif_minor_revision",
    ]
    predicate: Any = (
        [("entry_pdb_id", "==", normalized[0])]
        if len(normalized) == 1
        else [("entry_pdb_id", "in", normalized)]
    )
    rows = pd.read_parquet(path, columns=columns, filters=predicate)
    if rows["entry_pdb_id"].duplicated().any():
        duplicates = sorted(
            rows.loc[rows["entry_pdb_id"].duplicated(), "entry_pdb_id"].unique()
        )
        raise ValueError(f"duplicate source revisions in {path}: {duplicates}")
    revisions = {
        str(row.entry_pdb_id): (
            int(row.source_mmcif_major_revision),
            int(row.source_mmcif_minor_revision),
        )
        for row in rows.itertuples(index=False)
    }
    missing = sorted(set(normalized).difference(revisions))
    if missing:
        raise KeyError(f"source mmCIF revisions missing from {path}: {missing}")
    return revisions


def get_pdb_revision(
    pdb_or_system_id: str,
    *,
    manifest_path: Path | str | None = None,
) -> PDBRevision:
    """Return the exact source mmCIF revision pinned by this release."""
    pdb_id = _normalize_pdb_id(pdb_or_system_id)
    return _get_pdb_revisions([pdb_id], manifest_path=manifest_path)[pdb_id]


def _extended_pdb_id(pdb_id: str) -> str:
    """Return the wwPDB extended accession used by the versioned archive."""
    return pdb_id if pdb_id.startswith("pdb_") else f"pdb_0000{pdb_id}"


def pdb_mmcif_cache_path(
    pdb_or_system_id: str,
    *,
    cache_dir: Path | str | None = None,
    revision: PDBRevision | None = None,
    manifest_path: Path | str | None = None,
) -> Path:
    """Return the deterministic local path for a version-pinned PDB mmCIF."""
    pdb_id = _normalize_pdb_id(pdb_or_system_id)
    if revision is None:
        revision = get_pdb_revision(pdb_id, manifest_path=manifest_path)
    if cache_dir is None:
        cfg = get_config()
        cache_dir = Path(cfg.data.plinder_dir) / cfg.data.source_mmcifs
    major, minor = revision
    return Path(cache_dir) / pdb_id[-3:-1] / f"{pdb_id}_v{major}-{minor}.cif.gz"


def _mmcif_revision(path: Path) -> PDBRevision:
    """Read the latest structure-model revision from a compressed mmCIF."""
    try:
        with gzip.open(path, "rt", encoding="utf-8") as handle:
            cif_file = CIFFile.read(handle)
    except (
        OSError,
        EOFError,
        UnicodeDecodeError,
        DeserializationError,
        InvalidFileError,
    ) as exc:
        raise ValueError(f"invalid gzip-compressed mmCIF: {path}") from exc
    if not cif_file:
        raise ValueError(f"downloaded mmCIF contains no data block: {path}")
    block = next(iter(cif_file.values()))
    category_name = "pdbx_audit_revision_history"
    if category_name not in block:
        raise ValueError(f"mmCIF has no {category_name} category: {path}")
    category = block[category_name]
    if not {"major_revision", "minor_revision"}.issubset(category):
        raise ValueError(f"mmCIF has no revision columns: {path}")
    major = category["major_revision"].as_array()
    minor = category["minor_revision"].as_array()
    indices: Iterable[int] = range(len(major))
    if "data_content_type" in category:
        content_type = category["data_content_type"].as_array()
        indices = [
            index
            for index, value in enumerate(content_type)
            if str(value).lower() == "structure model"
        ]
    revisions = [(int(major[index]), int(minor[index])) for index in indices]
    if not revisions:
        raise ValueError(f"mmCIF contains no structure-model revision: {path}")
    return max(revisions)


def _validate_gzipped_mmcif(
    path: Path,
    *,
    expected_revision: PDBRevision | None = None,
) -> None:
    """Reject corrupt/non-mmCIF data and revision mismatches."""
    try:
        with gzip.open(path, "rb") as handle:
            prefix = handle.read(4096).lstrip()
    except (OSError, EOFError) as exc:
        raise ValueError(f"invalid gzip-compressed mmCIF: {path}") from exc
    if not prefix.startswith(b"data_"):
        raise ValueError(f"downloaded file is not an mmCIF: {path}")
    if expected_revision is not None:
        actual_revision = _mmcif_revision(path)
        if actual_revision != expected_revision:
            raise ValueError(
                f"source mmCIF revision mismatch for {path}: expected "
                f"{expected_revision[0]}.{expected_revision[1]}, found "
                f"{actual_revision[0]}.{actual_revision[1]}"
            )


def _download_versioned_mmcif(
    *,
    url: str,
    destination: Path,
    pdb_id: str,
    revision: PDBRevision,
) -> None:
    """Download one mmCIF atomically, retrying only transient failures."""
    major, minor = revision
    transient_errors = (
        requests.exceptions.Timeout,
        requests.exceptions.ConnectionError,
        requests.exceptions.ChunkedEncodingError,
        requests.exceptions.ContentDecodingError,
    )
    for attempt in range(1, SOURCE_MMCIF_DOWNLOAD_ATTEMPTS + 1):
        temporary: Path | None = None
        retry_reason: str | None = None
        try:
            response = requests.get(url, stream=True, timeout=(10, 120))
            if response.status_code == 404:
                raise FileNotFoundError(
                    f"wwPDB versioned archive has no source mmCIF for {pdb_id} "
                    f"revision {major}.{minor}: {url}"
                )
            if response.status_code in TRANSIENT_HTTP_STATUSES:
                if attempt == SOURCE_MMCIF_DOWNLOAD_ATTEMPTS:
                    response.raise_for_status()
                    raise RuntimeError(
                        f"HTTP {response.status_code} while downloading {url}"
                    )
                retry_reason = f"HTTP {response.status_code}"
            else:
                response.raise_for_status()
                with NamedTemporaryFile(
                    dir=destination.parent,
                    prefix=f".{destination.name}.",
                    suffix=".tmp",
                    delete=False,
                ) as handle:
                    temporary = Path(handle.name)
                    for chunk in response.iter_content(chunk_size=1024 * 1024):
                        if chunk:
                            handle.write(chunk)
                _validate_gzipped_mmcif(
                    temporary,
                    expected_revision=revision,
                )
                temporary.replace(destination)
                return
        except transient_errors as exc:
            if attempt == SOURCE_MMCIF_DOWNLOAD_ATTEMPTS:
                raise
            retry_reason = repr(exc)
        finally:
            if temporary is not None and temporary.exists():
                temporary.unlink()

        wait_seconds = 2 ** (attempt - 1)
        LOG.warning(
            f"transient source mmCIF download failure for {pdb_id} "
            f"({retry_reason}); retrying in {wait_seconds}s "
            f"[{attempt}/{SOURCE_MMCIF_DOWNLOAD_ATTEMPTS}]"
        )
        sleep(wait_seconds)

    raise RuntimeError(f"failed to download source mmCIF from {url}")


def get_pdb_mmcif(
    pdb_or_system_id: str,
    *,
    cache_dir: Path | str | None = None,
    force_update: bool = False,
    revision: PDBRevision | None = None,
    manifest_path: Path | str | None = None,
    base_url: str = WWPDB_VERSIONED_MMCIF_URL,
) -> Path:
    """Return the release-pinned source mmCIF, downloading it when allowed.

    Online calls populate the release-local cache atomically.  In offline mode
    the function performs no network or filesystem writes: a valid cached file
    is returned and a missing or corrupt file raises ``FileNotFoundError``.
    """
    pdb_id = _normalize_pdb_id(pdb_or_system_id)
    if revision is None:
        revision = get_pdb_revision(pdb_id, manifest_path=manifest_path)
    destination = pdb_mmcif_cache_path(
        pdb_id,
        cache_dir=cache_dir,
        revision=revision,
    )
    cached = False
    if destination.is_file():
        try:
            _validate_gzipped_mmcif(
                destination,
                expected_revision=revision,
            )
            cached = True
        except ValueError:
            cached = False
    if cached and (not force_update or is_offline()):
        return destination
    if is_offline():
        raise FileNotFoundError(
            f"No valid cached source mmCIF for {pdb_id} at {destination}. "
            "Disable PLINDER_OFFLINE/PLINDER_OFFLINE_MODE and call "
            "get_pdb_mmcif() or download_pdb_mmcifs() on an online node first."
        )

    destination.parent.mkdir(parents=True, exist_ok=True)
    major, minor = revision
    extended_id = _extended_pdb_id(pdb_id)
    filename = f"{extended_id}_xyz_v{major}-{minor}.cif.gz"
    url = f"{base_url.rstrip('/')}/{pdb_id[-3:-1]}/{extended_id}/{filename}"
    _download_versioned_mmcif(
        url=url,
        destination=destination,
        pdb_id=pdb_id,
        revision=revision,
    )
    LOG.info(
        f"cached wwPDB source mmCIF {pdb_id} revision {major}.{minor} at {destination}"
    )
    return destination


def download_pdb_mmcifs(
    pdb_or_system_ids: Iterable[str],
    *,
    cache_dir: Path | str | None = None,
    force_update: bool = False,
    max_workers: int = 8,
    manifest_path: Path | str | None = None,
) -> dict[str, Path]:
    """Prefetch a release-pinned subset of PDB mmCIFs for offline use."""
    if max_workers < 1:
        raise ValueError("max_workers must be at least 1")
    pdb_ids = list(
        dict.fromkeys(_normalize_pdb_id(value) for value in pdb_or_system_ids)
    )
    revisions = _get_pdb_revisions(pdb_ids, manifest_path=manifest_path)

    def download(pdb_id: str) -> Path:
        return get_pdb_mmcif(
            pdb_id,
            cache_dir=cache_dir,
            force_update=force_update,
            revision=revisions[pdb_id],
        )

    with ThreadPoolExecutor(max_workers=min(max_workers, len(pdb_ids) or 1)) as pool:
        downloaded = pool.map(download, pdb_ids)
    return dict(zip(pdb_ids, downloaded))


def retry(func: Callable[..., T]) -> Callable[..., T]:
    @wraps(func)
    def inner(*args: Any, **kwargs: Any) -> T:
        name = func.__name__
        mod = func.__module__
        log = setup_logger(".".join([mod, name]))
        retries = 5
        exc = None
        for i in range(1, retries + 1):
            try:
                return func(*args, **kwargs)
            except Exception as e:
                wait = 2**i
                log.error(f"failed: {repr(e)}, retry in: {wait}s")
                exc = e
                sleep(wait)
        raise Exception(f"Timeout error {exc}")

    return inner


@retry
def download_alphafold_cif_file(
    uniprot_id: str,
    output_folder: Path,
    url: str = "https://alphafold.ebi.ac.uk/files",
    force_update: bool = False,
) -> Optional[Path]:
    cif_file_path = output_folder / f"AF-{uniprot_id}-F1-model_v4.cif"
    if not cif_file_path.is_file() or force_update:
        resp = requests.get(f"{url}/{cif_file_path.name}")
        if resp.status_code == 404:
            LOG.info(f"UniProt ID {uniprot_id} not in AlphaFold database")
            return None
        resp.raise_for_status()
        with open(cif_file_path, "w") as f:
            f.write(resp.text)
    return cif_file_path


@retry
def download_pdb_chain_cif_file(pdb_id: str, chain_id: str, filename: Path) -> Path:
    structure = get_structure(
        CIFFile.read(
            fetch(
                pdb_ids=pdb_id,
                format="cif",
                overwrite=False,
            )
        ),
        model=1,
        use_author_fields=False,
        include_bonds=True,
    )
    write_file = CIFFile()
    set_structure(write_file, structure[structure.chain_id == chain_id])
    write_file.write(filename.as_posix())
    return filename
