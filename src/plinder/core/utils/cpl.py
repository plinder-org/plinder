"""Public dataset access with cloudpathlib and a local immutable-release cache."""

import os
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import Callable, Iterable, TypeVar
from urllib.parse import urlsplit

from tqdm import tqdm

from plinder.core.utils.config import get_config
from plinder.core.utils.r2 import ReleaseClient, checked_key, manifest

T = TypeVar("T")


def thread_pool(func: Callable[[T], None], items: Iterable[T]) -> None:
    with ThreadPoolExecutor(max_workers=8) as pool:
        list(pool.map(func, items))


def _get_client() -> ReleaseClient:
    cfg = get_config().data
    if (cfg.plinder_bucket, cfg.plinder_release, cfg.plinder_iteration) != (
        "plinder",
        "2024-06",
        "v2",
    ):
        raise ValueError("R2 downloads support only plinder/2024-06/v2")
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
        raise ValueError("PLINDER_MIRROR_URL must be an HTTP(S) bucket URL")
    if cfg.force_update:
        manifest.cache_clear()
    return ReleaseClient(remote)


def _download(client: ReleaseClient, paths: list[Path], force_progress: bool) -> None:
    cfg = get_config().data
    root = Path(cfg.plinder_dir).resolve()

    def fetch(path: Path) -> None:
        key = checked_key(path.resolve().relative_to(root).as_posix())
        if key not in client.records:
            raise FileNotFoundError(f"File absent from release manifest: {key}")
        if (
            cfg.force_update
            or not path.is_file()
            or path.stat().st_size != client.records[key]["size"]
        ):
            client.path(key).download_to(path)
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
    if paths and not os.getenv("PLINDER_OFFLINE"):
        _download(_get_client(), paths, force_progress)


def get_plinder_path(
    *, rel: str = "", download: bool = True, force_progress: bool = False
) -> Path:
    root = Path(get_config().data.plinder_dir)
    local = root / checked_key(rel)
    if os.getenv("PLINDER_OFFLINE") or not download:
        return local
    client = _get_client()
    remote = client.path(rel)
    if not remote.exists():
        raise FileNotFoundError(f"Path absent from release manifest: {rel}")
    files = (
        [remote] if remote.is_file() else [p for p in remote.rglob("*") if p.is_file()]
    )
    _download(client, [root / client.key(p) for p in files], force_progress)
    return local


def list_zip_paths(rel: str) -> list[Path]:
    root = Path(get_config().data.plinder_dir)
    rel = checked_key(rel)
    if os.getenv("PLINDER_OFFLINE"):
        return sorted((root / rel).glob("*.zip"))
    client = _get_client()
    return [root / client.key(p) for p in client.path(rel).glob("*.zip")]


def get_plinder_paths(*, paths: list[Path]) -> list[Path]:
    download_paths(paths=paths)
    return paths
