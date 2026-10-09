"""Public dataset access with cloudpathlib and a local release cache."""

import os
import re
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import Callable, Iterable, TypeVar
from urllib.parse import urlsplit

from tqdm import tqdm

from plinder.core.utils.config import get_config
from plinder.core.utils.http import ReleaseClient, checked_key

T = TypeVar("T")


def is_offline() -> bool:
    """Return whether remote access is disabled.

    ``PLINDER_OFFLINE`` remains the canonical setting.  The
    ``PLINDER_OFFLINE_MODE`` alias is accepted for compatibility with cluster
    launch environments that use the longer name.
    """
    values = (
        os.getenv("PLINDER_OFFLINE"),
        os.getenv("PLINDER_OFFLINE_MODE"),
    )
    return any(
        value is not None
        and value.strip().lower() not in {"", "0", "false", "no", "off"}
        for value in values
    )


def thread_pool(func: Callable[[T], None], items: Iterable[T]) -> None:
    with ThreadPoolExecutor(max_workers=8) as pool:
        list(pool.map(func, items))


def _get_client() -> ReleaseClient:
    cfg = get_config().data
    if cfg.plinder_release_number or not re.fullmatch(
        r"\d{4}-(0[1-9]|1[0-2])", cfg.plinder_release
    ):
        raise ValueError(
            "HTTP downloads support only YYYY-MM releases without a release number"
        )
    remote = str(cfg.plinder_remote).rstrip("/")
    url = urlsplit(remote)
    if (
        url.scheme not in {"http", "https"}
        or not url.netloc
        or url.query
        or url.fragment
        or url.username
        or url.password
    ):
        raise ValueError("PLINDER_MIRROR_URL must be an HTTP(S) file-server URL")
    return ReleaseClient(remote)


def _download(client: ReleaseClient, paths: list[Path], force_progress: bool) -> None:
    cfg = get_config().data
    root = Path(cfg.plinder_dir).resolve()

    def fetch(path: Path) -> None:
        key = checked_key(path.resolve().relative_to(root).as_posix())
        remote = client.path(key)
        if not remote.exists():
            raise FileNotFoundError(f"File absent from release: {key}")
        metadata = client.metadata(remote)
        if (
            cfg.force_update
            or not path.is_file()
            or path.stat().st_size != metadata.size
            or metadata.modified is None
            or path.stat().st_mtime != metadata.modified
        ):
            remote.download_to(path)
            if path.suffix == ".zip":
                path.with_name(path.stem + "_done").unlink(missing_ok=True)

    with ThreadPoolExecutor(max_workers=8) as pool:
        for _ in tqdm(
            pool.map(fetch, paths),
            total=len(paths),
            unit="file",
            disable=not force_progress,
        ):
            pass


def download_paths(*, paths: list[Path], force_progress: bool = False) -> None:
    """
    Download pre-determined paths from the release mirror concurrently. This
    is useful when we want to process a pre-determined subset of the data
    rather than all of the contents of the dataset.

    Parameters
    ----------
    paths : list[Path]
        the local paths to download, resolved against the release cache root
    force_progress : bool, default=False
        if True, always display a progress bar
    """
    if paths and not is_offline():
        _download(_get_client(), paths, force_progress)


def _prune_unpublished(directory: Path, published: set[str]) -> None:
    """Remove cached Parquet files that the release no longer lists.

    Directory artifacts are read as whole Parquet datasets, so a file dropped
    from an updated release would otherwise be read alongside its replacement.
    """
    root = Path(get_config().data.plinder_dir)
    for path in directory.rglob("*.parquet"):
        if path.is_symlink() or not path.is_file():
            continue
        if path.relative_to(root).as_posix() not in published:
            path.unlink()


def get_plinder_path(
    *, rel: str = "", download: bool = True, force_progress: bool = False
) -> Path:
    """
    Get the local cache path for a file or directory in the plinder release,
    downloading it first if requested. This provides a convenient way to
    manage local file caching since it is automatically synced from the
    release mirror on access in case remote files change.

    Parameters
    ----------
    rel : str
        Relative path to the file or directory.
    download : bool, default=True
        if True, download the files
    force_progress : bool, default=False
        if True, force progress bar even if < 10 files

    Returns
    -------
    Path
        The local cache path.
    """
    root = Path(get_config().data.plinder_dir)
    local = root / checked_key(rel)
    if is_offline() or not download:
        return local
    local.resolve().relative_to(root.resolve())
    client = _get_client()
    remote = client.path(rel)
    if not remote.exists():
        raise FileNotFoundError(f"Path absent from release: {rel}")
    files = (
        [remote] if remote.is_file() else [p for p in remote.rglob("*") if p.is_file()]
    )
    keys = {client.key(p) for p in files}
    _download(client, [root / key for key in sorted(keys)], force_progress)
    if not remote.is_file():
        local.mkdir(parents=True, exist_ok=True)
        _prune_unpublished(local, keys)
    return local
