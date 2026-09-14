"""Manifest-backed listing and verified cloudpathlib downloads for public R2."""

import base64
import gzip
import hashlib
import json
from functools import lru_cache
from http.client import HTTPException
from os import PathLike
from pathlib import Path
from tempfile import TemporaryDirectory
from time import sleep
from typing import Iterator, TypedDict, TypeVar, cast
from urllib.error import HTTPError
from urllib.parse import quote, unquote

from cloudpathlib import CloudPath as BaseCloudPath
from cloudpathlib import HttpPath, HttpsClient, HttpsPath

CloudPathT = TypeVar("CloudPathT", bound=BaseCloudPath)


class ManifestRecord(TypedDict):
    key: str
    size: int
    md5: str


def checked_key(key: str) -> str:
    if key and ("\\" in key or any(p in {"", ".", ".."} for p in key.split("/"))):
        raise ValueError(f"Invalid dataset path: {key!r}")
    return key


@lru_cache(maxsize=4)
def manifest(remote: str) -> dict[str, ManifestRecord]:
    with TemporaryDirectory() as tmp:
        path = Path(tmp) / "manifest.jsonl.gz"
        path_class = HttpsPath if remote.startswith("https:") else HttpPath
        path_class(remote + "/manifest.jsonl.gz").download_to(path)
        data = path.read_bytes()
    records: dict[str, ManifestRecord] = {}
    for line in gzip.decompress(data).splitlines():
        row = json.loads(line)
        key = checked_key(row["key"])
        if not key or key in records or type(row["size"]) is not int or row["size"] < 0:
            raise ValueError("Invalid or duplicate manifest record")
        if len(base64.b64decode(row["md5"], validate=True)) != 16:
            raise ValueError("Invalid manifest MD5")
        records[key] = row
    return records


class ReleaseClient(HttpsClient):
    """Use the manifest where public HTTP lacks object listing; delegate transfers."""

    def __init__(self, remote: str):
        super().__init__()
        self.remote = remote
        self.records = manifest(remote)
        self.dir_matcher = lambda url: self.key(url) not in self.records

    def CloudPath(self, path: str | CloudPathT, *parts: str) -> CloudPathT:
        cls = HttpsPath if self.remote.startswith("https:") else HttpPath
        return cast(CloudPathT, cls(str(path), *parts, client=self))

    def path(self, key: str = "") -> HttpPath:
        return cast(
            HttpPath,
            self.CloudPath(
                self.remote + ("/" + quote(checked_key(key), safe="/") if key else "")
            ),
        )

    def key(self, path: str | BaseCloudPath) -> str:
        url = str(path)
        if url == self.remote:
            return ""
        if not url.startswith(self.remote + "/"):
            raise ValueError("Path is outside the release")
        return checked_key(unquote(url[len(self.remote) + 1 :]))

    def _exists(self, path: HttpPath) -> bool:
        key = self.key(path)
        return key in self.records or any(
            not key or k.startswith(key + "/") for k in self.records
        )

    def _list_dir(
        self, path: HttpPath, recursive: bool
    ) -> Iterator[tuple[HttpPath, bool]]:
        key = self.key(path)
        prefix = key + "/" if key else ""
        children = {}
        for name in self.records:
            if name.startswith(prefix):
                parts = name[len(prefix) :].split("/")
                for i in range(1, len(parts) + 1 if recursive else 2):
                    children[prefix + "/".join(parts[:i])] = i < len(parts)
        for name, is_dir in sorted(children.items()):
            yield self.path(name), is_dir

    def _download_file(self, path: HttpPath, local_path: str | PathLike[str]) -> Path:
        destination = Path(local_path)
        destination.parent.mkdir(parents=True, exist_ok=True)
        record = self.records[self.key(path)]
        with TemporaryDirectory(dir=destination.parent) as tmp:
            temporary = Path(tmp) / "download"
            for attempt in range(3):
                try:
                    super()._download_file(path, temporary)
                    if temporary.stat().st_size != record["size"]:
                        raise OSError("Incomplete dataset download")
                    break
                except (OSError, HTTPException) as exc:
                    if isinstance(exc, HTTPError) and exc.code not in {
                        408,
                        429,
                        500,
                        502,
                        503,
                        504,
                    }:
                        raise
                    if attempt == 2:
                        raise
                    sleep(2**attempt)
            digest = hashlib.md5()
            with temporary.open("rb") as stream:
                for chunk in iter(lambda: stream.read(1024 * 1024), b""):
                    digest.update(chunk)
            if base64.b64encode(digest.digest()).decode() != record["md5"]:
                raise ValueError("Dataset checksum mismatch")
            temporary.replace(destination)
        return destination
