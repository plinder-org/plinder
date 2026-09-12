#!/usr/bin/env python
# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
"""Prepare shared ingest reference data on a machine with internet access."""

from __future__ import annotations

import argparse
import os
import time
from pathlib import Path

from plinder.data.pipeline.tasks import download_alternative_datasets


def _is_truthy(value: str | None) -> bool:
    return value is not None and value.lower() in {"1", "true", "yes", "on"}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output_root", type=Path)
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()

    if _is_truthy(os.environ.get("PLINDER_OFFLINE")) or _is_truthy(
        os.environ.get("PLINDER_OFFLINE_MODE")
    ):
        raise RuntimeError(
            "reference preparation requires internet access; unset "
            "PLINDER_OFFLINE and PLINDER_OFFLINE_MODE"
        )

    started = time.perf_counter()
    download_alternative_datasets(
        data_dir=args.output_root.resolve(),
        threads=args.threads,
        force_update=args.force,
    )
    print(
        f"prepared reference data in {time.perf_counter() - started:.1f}s "
        f"under {args.output_root.resolve() / 'dbs'}"
    )


if __name__ == "__main__":
    main()
