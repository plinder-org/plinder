# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
import os
from concurrent.futures import ALL_COMPLETED, ThreadPoolExecutor, wait
from functools import wraps
from pathlib import Path
from time import sleep
from typing import Callable, Iterable, Optional, TypeVar, Union, overload

from cloudpathlib import CloudPath, GSClient, GSPath
from cloudpathlib.exceptions import OverwriteNewerLocalError
from google.cloud import storage
from omegaconf import DictConfig
from tqdm.contrib.concurrent import thread_map

from plinder.core.utils.config import get_config
from plinder.core.utils.dec import timeit
from plinder.core.utils.log import setup_logger

T = TypeVar("T")
LOG = setup_logger(__name__)
_CLIENTS: dict[str, GSClient] = {}


def _retry_decorator(retries: int) -> Callable[[Callable[..., T]], Callable[..., T]]:
    def decorator(func: Callable[..., T]) -> Callable[..., T]:
        @wraps(func)
        def inner(path: GSPath) -> T:
            exc: Optional[Exception] = None
            for i in range(1, retries + 1):
                try:
                    return func(path)
                except Exception as e:
                    LOG.error(f"retry {path}: {repr(e)}")
                    exc = e
                    sleep(2**i)
            raise Exception(f"Retry error {path}: {repr(exc)}") from exc

        return inner

    return decorator


@overload
def retry(f: Callable[..., T]) -> Callable[..., T]:
    ...


@overload
def retry(*, retries: int = 5) -> Callable[[Callable[..., T]], Callable[..., T]]:
    ...


def retry(
    f: Optional[Callable[..., T]] = None, retries: int = 5
) -> Union[Callable[..., T], Callable[[Callable[..., T]], Callable[..., T]]]:
    if f is None:
        return _retry_decorator(retries)
    return _retry_decorator(retries)(f)


def thread_pool(func: Callable[..., None], iter: Iterable[T]) -> None:
    with ThreadPoolExecutor() as executor:
        futures = [executor.submit(func, item) for item in iter]
        wait(futures, return_when=ALL_COMPLETED)
        for future in futures:
            exc = future.exception()
            if exc is not None:
                raise exc


@retry
def _quiet_ping(path: GSPath) -> None:
    if isinstance(path, CloudPath):
        LOG.debug(f"_ping: type(path)={path.__class__.__name__} local={path._local}")
        if not os.getenv("PLINDER_OFFLINE"):
            try:
                path.fspath
            except OverwriteNewerLocalError:
                if path._local.exists():
                    pass
                else:
                    raise
    else:
        LOG.debug(f"_ping: type(path)={type(path)} local={path}")


@timeit
def download_paths(*, paths: list[GSPath], force_progress: bool = False) -> None:
    """
    Download pre-determined paths from GCS concurrently. This is useful
    when we want to process a pre-determined subset of the data rather
    than all of the contents of the dataset.

    Parameters
    ----------
    paths : list[GSPath]
        the paths to download
    """
    if os.getenv("PLINDER_OFFLINE"):
        return
    if len(paths) > 10 or force_progress:
        thread_map(_quiet_ping, paths)
    else:
        thread_pool(_quiet_ping, paths)


def _get_client(cfg: DictConfig) -> GSClient:
    bucket = cfg.data.plinder_bucket
    client = _CLIENTS.get(bucket)
    if client is None:
        mount = cfg.data.plinder_mount
        if mount in ["/plinder", "/", ""]:
            mount = "/"
        if bucket == "plinder":
            client = GSClient(
                local_cache_dir=mount,
                storage_client=storage.Client.create_anonymous_client(),
            )
        else:
            client = GSClient(local_cache_dir=mount)
        _CLIENTS[bucket] = client
    return _CLIENTS[bucket]


def get_plinder_path(
    *, rel: str = "", download: bool = True, force_progress: bool = False
) -> Path:
    """
    Get a cloudpathlib path to a file or directory in the plinder bucket.
    This provides a convenient way to manage local file caching since it
    is automatically synced from the cloud on access in case remote files
    change.

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
    GSPath
        The cloudpathlib path.
    """
    cfg = get_config()
    client = _get_client(cfg)
    remote = cfg.data.plinder_remote
    if rel:
        remote += f"/{rel}"
    path = GSPath(remote, client=client)
    LOG.debug(f"get_plinder_path: remote={path}")
    if os.getenv("PLINDER_OFFLINE"):
        return Path(path._local)

    if download:
        paths = [path] if path.is_file() else list(path.rglob("*"))
        download_paths(paths=paths, force_progress=force_progress)
    try:
        return Path(path.fspath)
    except OverwriteNewerLocalError:
        if path._local.exists():
            return Path(path._local)
        else:
            raise


def get_plinder_paths(*, paths: list[Path]) -> list[Path]:
    cfg = get_config()
    client = _get_client(cfg)
    remote = GSPath(cfg.data.plinder_remote, client=client)
    LOG.debug(f"get_plinder_paths: remote={remote} npaths={len(paths)}")
    anypaths = [remote / path.relative_to(cfg.data.plinder_dir) for path in paths]
    if not os.getenv("PLINDER_OFFLINE"):
        download_paths(paths=anypaths)
    return anypaths
