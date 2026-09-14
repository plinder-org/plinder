"""Build the release's HTTP listing from audited GCS metadata, without cloud I/O."""

import argparse
import base64
import gzip
import hashlib
import json
from pathlib import Path


def build(inventory: Path, destination: Path) -> int:
    shards = json.loads((inventory / "manifest-files.json").read_text())
    records = {}
    for shard in shards:
        path = (inventory / shard["path"]).resolve()
        if not path.is_relative_to(inventory.resolve()):
            raise ValueError("Invalid inventory path")
        if hashlib.sha256(path.read_bytes()).hexdigest() != shard["sha256"]:
            raise ValueError("Inventory checksum mismatch")
        with gzip.open(path, "rt") as stream:
            for line in stream:
                r = json.loads(line)
                if not r["name"].startswith("2024-06/v2/"):
                    raise ValueError("Unexpected release")
                key = r["name"][len("2024-06/v2/") :]
                if key.endswith("/"):
                    continue
                if (
                    not key
                    or "\\" in key
                    or any(p in {"", ".", ".."} for p in key.split("/"))
                    or key in records
                ):
                    raise ValueError("Invalid or duplicate object key")
                size = int(r["size"])
                if size < 0 or len(base64.b64decode(r["md5Hash"], validate=True)) != 16:
                    raise ValueError("Invalid object metadata")
                records[key] = {"key": key, "size": size, "md5": r["md5Hash"]}
    destination.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(destination, "wt") as stream:
        for key in sorted(records):
            stream.write(json.dumps(records[key], separators=(",", ":")) + "\n")
    return len(records)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("inventory", type=Path)
    parser.add_argument("destination", type=Path)
    args = parser.parse_args()
    print(build(args.inventory, args.destination))
