# Release downloads from the Cameo file server

Public dataset reads for `YYYY-MM` releases (currently `2026-09`) use
[the Cameo file server](https://cameo3d.org/plinder/PLINDER-2026-09/) directly.
No S3 keys, Google credentials, or separately published manifest are required.
Releases with a release number, such as the legacy `2024-06/v2`, are rejected;
use an earlier `plinder` version for them.

`PLINDER_MIRROR_URL` may select another HTTP(S) file-server root containing
`PLINDER-<YYYY-MM>/` directories with nginx or Python-style directory indexes.
For example, setting it to `https://mirror.example/datasets` selects
`https://mirror.example/datasets/PLINDER-2026-09/` for the default release.
The local cache remains `<PLINDER_MOUNT>/<PLINDER_BUCKET>/<YYYY-MM>/`.

## Download and cache behavior

The consumer uses cloudpathlib HTTP paths. Directory indexes provide file
listings; HEAD responses provide exact sizes and Last-Modified times. Listings
and metadata are fetched afresh for each API download call, so additions,
removals, and hotfixes are visible without restarting Python. Directory links
are restricted to immediate children within the release.

Transfers stream into temporary files, validate HTTP Content-Length, and
atomically replace the destination. Interrupted transfers retry from the start,
with up to three attempts. They do not resume at a byte offset. Failed transfers
preserve an existing destination. Metadata and directory requests use the same
retry policy. Parallel download batches use eight workers and propagate errors.

The file server does not publish per-file checksums. Content-Length validation
detects truncated transfers; it does not detect same-size content corruption.

Downloaded files retain the server's Last-Modified time. Cached files are reused
when both size and modification time match. A changed server timestamp triggers
a new download even for a same-size hotfix. Files without Last-Modified are
fetched again on access. Publishers must update Last-Modified when changing a
file; a same-size edit preserving that timestamp requires a forced refresh
(`data.force_update=true`) or deletion of the affected cached file.
Refreshing a ZIP clears its extraction marker.

When a published release is updated in place, a directory download also removes
cached Parquet files under that directory that the server no longer lists.
Pruning happens only after the complete listing and all downloads succeed, so
readers do not mix removed fragments with their replacements.

`PLINDER_OFFLINE=true` bypasses all network access; unset it to return online.
`PLINDER_OFFLINE_MODE` is also accepted. Empty values and `0`, `false`, `no`,
and `off` leave networking enabled. Google libraries are optional for consumers
and available through `plinder[data]` for data generation.

## Focused validation

```bash
python -m pytest --noconftest tests/core/test_http.py tests/core/test_core_config.py -q -o addopts=''
```

The local HTTP-server tests exercise directory and direct-file access without a
manifest, encoded and nested paths, cache reuse, same-size server hotfixes,
forced refresh, removed Parquet fragments, offline behavior, path containment,
atomic replacement, and recovery from interrupted transfers and HTTP errors.
`--noconftest` excludes the repository-wide fixtures; full project CI also runs
the release APIs and executes all documentation notebooks against the default
public server.
