"""Directory listing and atomic downloads from the public PLINDER file server."""

import os
import shutil
from dataclasses import dataclass
from email.utils import parsedate_to_datetime
from html.parser import HTMLParser
from http.client import HTTPException
from os import PathLike
from pathlib import Path
from tempfile import TemporaryDirectory
from time import sleep
from typing import Callable, Iterator, TypeVar, cast
from urllib.error import HTTPError
from urllib.parse import quote, unquote, urlsplit
from urllib.request import BaseHandler, Request

from cloudpathlib import CloudPath as BaseCloudPath
from cloudpathlib import HttpPath, HttpsClient, HttpsPath

CloudPathT = TypeVar("CloudPathT", bound=BaseCloudPath)
T = TypeVar("T")


def checked_key(key: str) -> str:
    if key and ("\\" in key or any(p in {"", ".", ".."} for p in key.split("/"))):
        raise ValueError(f"Invalid dataset path: {key!r}")
    return key


class _TimeoutHandler(BaseHandler):
    def http_request(self, request: Request) -> Request:
        # Bound connection/read inactivity without limiting total transfer time.
        request.timeout = 60
        return request

    https_request = http_request


def _retry(operation: Callable[[], T]) -> T:
    for attempt in range(3):
        try:
            return operation()
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
    raise AssertionError("unreachable")


class _DirectoryIndex(HTMLParser):
    def __init__(self) -> None:
        super().__init__()
        self.links: list[str] = []
        self.title = ""
        self.in_title = False

    def handle_starttag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
        if tag == "title":
            self.in_title = True
        if tag == "a":
            self.links.extend(
                value for name, value in attrs if name == "href" and value
            )

    def handle_endtag(self, tag: str) -> None:
        if tag == "title":
            self.in_title = False

    def handle_data(self, data: str) -> None:
        if self.in_title:
            self.title += data


@dataclass(frozen=True)
class FileMetadata:
    size: int
    modified: float | None
    is_dir: bool


class ReleaseClient(HttpsClient):
    """Read a fresh directory index and HTTP metadata for each cache access."""

    def __init__(self, remote: str):
        super().__init__(auth=_TimeoutHandler())
        self.remote = remote.rstrip("/")
        # These caches live only for one get_plinder_path/download_paths call.
        self._metadata: dict[str, FileMetadata] = {}
        self._listed: dict[str, bool] = {}
        self.dir_matcher = self._is_directory

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
        url = str(path).removesuffix("/")
        if url == self.remote:
            return ""
        if not url.startswith(self.remote + "/"):
            raise ValueError("Path is outside the release")
        return checked_key(unquote(url[len(self.remote) + 1 :]))

    def metadata(self, path: HttpPath) -> FileMetadata:
        key = self.key(path)
        if key not in self._metadata:

            def read() -> FileMetadata:
                request = Request(
                    str(path), method="HEAD", headers={"Cache-Control": "no-cache"}
                )
                with self.opener.open(request) as response:
                    if self.key(response.geturl()) != key:
                        raise ValueError(
                            "File server redirected outside the requested path"
                        )
                    is_dir = response.geturl().endswith("/")
                    size = response.headers.get("Content-Length")
                    if not is_dir and size is None:
                        raise ValueError("File server must supply Content-Length")
                    modified = response.headers.get("Last-Modified")
                    return FileMetadata(
                        int(size or 0),
                        parsedate_to_datetime(modified).timestamp()
                        if modified
                        else None,
                        is_dir,
                    )

            self._metadata[key] = _retry(read)
        return self._metadata[key]

    def _is_directory(self, url: str) -> bool:
        key = self.key(url)
        if key in self._listed:
            return self._listed[key]
        return self.metadata(self.path(key)).is_dir

    def _exists(self, path: HttpPath) -> bool:
        if self.key(path) in self._listed:
            return True
        try:
            self.metadata(path)
        except HTTPError as exc:
            if exc.code == 404:
                return False
            raise
        return True

    def _list_dir(
        self, path: HttpPath, recursive: bool
    ) -> Iterator[tuple[HttpPath, bool]]:
        key = self.key(path)
        url = str(self.path(key)) + "/"

        def read() -> str:
            request = Request(url, headers={"Cache-Control": "no-cache"})
            with self.opener.open(request) as response:
                if self.key(response.geturl()) != key:
                    raise ValueError(
                        "File server redirected outside the requested directory"
                    )
                return cast(bytes, response.read()).decode("utf-8")

        index = _DirectoryIndex()
        index.feed(_retry(read))
        if not index.title.startswith(("Index of ", "Directory listing for ")):
            raise ValueError("Expected an HTTP directory index")
        children: dict[str, bool] = {}
        for href in index.links:
            link = urlsplit(href)
            if link.scheme or link.netloc or link.query or link.fragment:
                continue
            name = unquote(link.path.removesuffix("/"))
            # Only immediate relative children: ignore parent/sort/navigation links.
            if not name or name in {".", ".."} or "/" in name or "\\" in name:
                continue
            child = checked_key(f"{key}/{name}" if key else name)
            children[child] = link.path.endswith("/")
        self._listed.update(children)
        for child, is_dir in sorted(children.items()):
            child_path = self.path(child)
            yield child_path, is_dir
            if recursive and is_dir:
                yield from self._list_dir(child_path, recursive=True)

    def _download_file(self, path: HttpPath, local_path: str | PathLike[str]) -> Path:
        destination = Path(local_path)
        destination.parent.mkdir(parents=True, exist_ok=True)
        with TemporaryDirectory(dir=destination.parent) as tmp:
            temporary = Path(tmp) / "download"

            def download() -> None:
                request = Request(str(path), headers={"Cache-Control": "no-cache"})
                with self.opener.open(request) as response:
                    if self.key(response.geturl()) != self.key(path):
                        raise ValueError(
                            "File server redirected outside the requested path"
                        )
                    size = response.headers.get("Content-Length")
                    if size is None:
                        raise ValueError("File server must supply Content-Length")
                    with temporary.open("wb") as stream:
                        shutil.copyfileobj(response, stream)
                    if temporary.stat().st_size != int(size):
                        raise OSError("Incomplete dataset download")
                    modified = response.headers.get("Last-Modified")
                    if modified:
                        timestamp = parsedate_to_datetime(modified).timestamp()
                        os.utime(temporary, (timestamp, timestamp))

            _retry(download)
            temporary.replace(destination)
        return destination
