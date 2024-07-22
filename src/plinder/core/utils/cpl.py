# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
from functools import wraps
from time import sleep
from typing import Callable, Optional, TypeVar, Union, overload

from cloudpathlib import AnyPath, GSClient
from tqdm.contrib.concurrent import thread_map

from plinder.core.utils.config import get_config
from plinder.core.utils.dec import timeit
from plinder.core.utils.log import setup_logger

T = TypeVar("T")
LOG = setup_logger(__name__)


def _retry_decorator(retries: int) -> Callable[[Callable[..., T]], Callable[..., T]]:
    def decorator(func: Callable[..., T]) -> Callable[..., T]:
        @wraps(func)
        def inner(path: AnyPath) -> T:
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


@retry
@timeit
def ping(path: AnyPath) -> None:
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

    @retry
    def _ping(path: AnyPath) -> None:
        path.fspath

    thread_map(_ping, list(root.rglob("*")))


def get_plinder_path(*, rel: str) -> AnyPath:
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
    AnyPath
        The cloudpathlib path.
    """
    cfg = get_config()
    root = cfg.data.plinder_mount
    if hasattr(cfg, "ingest"):
        root = cfg.ingest.plinder_mount
    if root == "/plinder":
        root = "/"
    elif root == "":
        root = "/"
    client = GSClient(local_cache_dir=root)
    remote = f"{cfg.data.plinder_remote}/{rel}"
    path = AnyPath(remote, client=client)
    LOG.debug(f"gspath: remote={remote}")
    ping(path)
    LOG.debug(f"fspath: local={path.fspath}")
    return path
