# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
from concurrent.futures import ThreadPoolExecutor
from functools import wraps
from pathlib import Path
from time import sleep
from typing import Any, Callable, Optional, TypeVar

from google.cloud.storage.bucket import Bucket
from google.cloud.storage.client import Client
from omegaconf import DictConfig
from tqdm import tqdm

from plinder.core.utils.dec import timeit
from plinder.core.utils.log import setup_logger

BUCKET = "plinder"
BUCKETS: dict[str, Bucket] = {}
LOG = setup_logger(__name__)
T = TypeVar("T")


def retry(f: Optional[Callable[..., T]] = None, retries: int = 5) -> Callable[..., T]:
    def decorator(func: Callable[..., T]) -> Callable[..., T]:
        @wraps(func)
        def inner(**kwargs: Any) -> T:
            bucket = kwargs.get("bucket", None)
            if bucket is None:
                cfg = kwargs.get("cfg", None)
                if cfg is None:
                    name = kwargs.get("bucket_name")
                    if name is not None:
                        client = Client()
                        bucket = client.bucket(name)
                        BUCKETS[name] = bucket
                    else:
                        client = Client.create_anonymous_client()
                        bucket = client.bucket(BUCKET)
                        BUCKETS[BUCKET] = bucket
                else:
                    if cfg.data.plinder_bucket in BUCKETS:
                        bucket = BUCKETS[cfg.data.plinder_bucket]
                    else:
                        if cfg.data.plinder_bucket == BUCKET:
                            client = Client.create_anonymous_client()
                        else:
                            client = Client()
                        bucket = client.bucket(cfg.data.plinder_bucket)
                        BUCKETS[cfg.data.plinder_bucket] = bucket
            exc = None
            for i in range(1, retries + 1):
                try:
                    return func(bucket=bucket, **kwargs)
                except Exception as e:
                    LOG.error(f"retry: {e}")
                    exc = e
                    sleep(2**i)
            raise Exception(f"Timeout error {exc}")

        return inner

    return decorator(f) if f else decorator  # type: ignore


@retry
def download_as_str(
    *,
    gcs_path: str,
    cfg: Optional[DictConfig] = None,
    bucket_name: Optional[str] = None,
    bucket: Optional["Bucket"] = None,
) -> str:
    assert bucket is not None
    blob = bucket.blob("/".join(Path(gcs_path).parts[2:]))
    return str(blob.download_as_string().decode("utf8"))


@retry
def download_to_file(
    *,
    gcs_path: str,
    local_path: str,
    cfg: Optional[DictConfig] = None,
    bucket_name: Optional[str] = None,
    bucket: Optional["Bucket"] = None,
) -> None:
    assert bucket is not None
    blob = bucket.blob("/".join(Path(gcs_path).parts[2:]))
    blob.download_to_filename(local_path)


@retry
def list_dir(
    *,
    gcs_path: str,
    cfg: Optional[DictConfig] = None,
    bucket_name: Optional[str] = None,
    bucket: Optional["Bucket"] = None,
) -> list[str]:
    assert bucket is not None
    prefix = "/".join(Path(gcs_path.rstrip("/")).parts[2:])
    return [
        f"gs://{bucket.name}/{blob.name}"
        for blob in bucket.list_blobs(prefix=f"{prefix}/")
    ]


@retry
def upload_from_file(
    *,
    local_path: str,
    gcs_path: str,
    cfg: Optional[DictConfig] = None,
    bucket_name: Optional[str] = None,
    bucket: Optional["Bucket"] = None,
) -> None:
    assert bucket is not None
    blob = bucket.blob("/".join(Path(gcs_path).parts[2:]))
    with open(local_path, "rb") as f:
        blob.upload_from_file(f)


@timeit
def download_many(
    *,
    gcs_paths: list[str],
    local_paths: list[str],
    cfg: Optional[DictConfig] = None,
    bucket_name: Optional[str] = None,
) -> None:
    if not len(gcs_paths):
        LOG.info("No files to download")
        return
    if not len(gcs_paths) == len(local_paths):
        raise ValueError("gcs_paths and local_paths must have the same length")

    def ordered(gcs_path: str, local_path: str) -> None:
        download_to_file(
            gcs_path=gcs_path, local_path=local_path, cfg=cfg, bucket_name=bucket_name
        )

    with ThreadPoolExecutor() as executor:
        res = executor.map(ordered, gcs_paths, local_paths)
        if len(gcs_paths) > 20:
            list(tqdm(res, total=len(gcs_paths)))
        else:
            LOG.info(f"downloading {len(gcs_paths)} files")


def download_if_not_exists(
    *,
    local_paths: list[Path],
    remote_root: Path,
    cfg: Optional[DictConfig] = None,
    bucket_name: Optional[str] = None,
) -> None:
    gcs_paths = []
    need_download = []
    for local_path in local_paths:
        if not local_path.is_file():
            need_download.append(local_path.as_posix())
            gcs_paths.append((remote_root / local_path.name).as_posix())
    if len(need_download):
        download_many(
            gcs_paths=gcs_paths,
            local_paths=need_download,
            cfg=cfg,
        )
