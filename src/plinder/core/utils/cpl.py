# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
from concurrent.futures import ALL_COMPLETED, ThreadPoolExecutor, wait
from functools import wraps
from pathlib import Path
from time import sleep
from typing import Callable, Iterable, Optional, TypeVar, Union, overload

from cloudpathlib import CloudPath, GSClient, GSPath
from omegaconf import DictConfig
from tqdm.contrib.concurrent import thread_map

from plinder.core.utils.config import get_config
from plinder.core.utils.dec import timeit
from plinder.core.utils.log import setup_logger

T = TypeVar("T")
LOG = setup_logger(__name__)


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


def thread_pool(func: Callable[..., T], iter: Iterable[T]) -> None:
    with ThreadPoolExecutor() as executor:
        futures = [executor.submit(func, item) for item in iter]
        wait(futures, return_when=ALL_COMPLETED)
        for future in futures:
            exc = future.exception()
            if exc is not None:
                raise exc


@retry
def _quiet_ping(path: GSPath) -> None:
    LOG.debug(f"_ping: type(path)={type(path)} local={path.fspath}")
    if isinstance(path, CloudPath):
        path.fspath


@timeit
def download_many(*, rel: str) -> None:
    """
    Download many files from GCS concurrently. This is useful when
    it is known that we need access to many files in a single operation.
    If the files were already downloaded, you still pay the cost of
    checking if the remote file has drifted from the local file even
    if they are identical.

    Parameters
    ----------
    rel : str
        Relative path to the files to download.
    """
    root = get_plinder_path(rel=rel)

    paths = list(root.rglob("*"))

    if len(paths) > 10:
        thread_map(_quiet_ping, paths)
    else:
        thread_pool(_quiet_ping, paths)


@timeit
def download_paths(*, paths: list[GSPath]) -> None:
    """
    Download pre-determined paths from GCS concurrently. This is useful
    when we want to process a pre-determined subset of the data rather
    than all of the contents of the dataset.

    Parameters
    ----------
    paths : list[GSPath]
        the paths to download
    """
    if len(paths) > 10:
        thread_map(_quiet_ping, paths)
    else:
        thread_pool(_quiet_ping, paths)


def _get_fsroot(cfg: DictConfig) -> str:
    """
    Kludge mainly for dealing with the NFS
    """
    root = cfg.data.plinder_mount
    if root in ["/plinder", "/", ""]:
        root = "/"
    return root


def get_plinder_path(*, rel: str = "") -> CloudPath:
    """
    Get a cloudpathlib path to a file or directory in the plinder bucket.
    This provides a convenient way to manage local file caching since it
    is automatically synced from the cloud on access in case remote files
    change.

    Parameters
    ----------
    rel : str
        Relative path to the file or directory.

    Returns
    -------
    GSPath
        The cloudpathlib path.
    """
    cfg = get_config()
    root = _get_fsroot(cfg)
    client = GSClient(local_cache_dir=root)
    remote = cfg.data.plinder_remote
    if rel:
        remote = f"{cfg.data.plinder_remote}/{rel}"
    path = GSPath(remote, client=client)
    LOG.debug(f"get_plinder_path: remote={path} root={root}")
    return path


def get_plinder_paths(*, paths: list[Path]) -> list[CloudPath]:
    cfg = get_config()
    root = _get_fsroot(cfg)
    client = GSClient(local_cache_dir=root)
    remote = GSPath(cfg.data.plinder_remote, client=client)
    LOG.debug(f"get_plinder_paths: remote={remote} root={root}")
    anypaths = [remote / path.relative_to(root) for path in paths]
    return anypaths
