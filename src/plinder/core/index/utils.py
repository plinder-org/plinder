# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
from __future__ import annotations

from argparse import ArgumentParser
from textwrap import dedent
from time import time

from plinder.core.release import RELEASE_PATHS, RELEASE_TABLES
from plinder.core.utils import cpl
from plinder.core.utils.config import get_config
from plinder.core.utils.log import setup_logger

LOG = setup_logger(__name__)


_DOWNLOAD_GROUPS = (
    (
        "index tables",
        tuple(
            dict.fromkeys(str(table["artifact"]) for table in RELEASE_TABLES.values())
        ),
        False,
    ),
    (
        "representative-cover tables",
        ("ligand_sampling", "interface_sampling"),
        False,
    ),
    (
        "complete similarity exports",
        ("ligand_similarity_scores", "interface_similarity_scores"),
        True,
    ),
    ("canonical ligand archives", ("ligand_archives",), True),
    ("ligand similarities", ("ligand_scores",), True),
    ("protein-interface similarities", ("interface_scores",), True),
    ("mapped protein alignments", ("alignments",), True),
    ("custom-scoring search databases", ("search_databases",), True),
)


def download_plinder_cmd(args: list[str] | None = None) -> None:
    """
    Download published PLINDER release artifacts for the current configuration.

    Source PDB mmCIFs are not included. They are fetched and cached per PDB
    entry when reconstruction needs them.

    Note that even though this is wrapped in a progress bar, the estimated
    completion time can vary wildly as it iterates over larger files vs.
    smaller ones.
    """
    t0 = time()
    parser = ArgumentParser(usage=download_plinder_cmd.__doc__)
    parser.add_argument(
        "--release",
        default=None,
        help="ingest month for the PLINDER release (YYYY-MM)",
    )
    parser.add_argument(
        "--release-number",
        default=None,
        help="numbered release within the ingest month",
    )
    parser.add_argument("-y", "--yes", action="store_true", help="skip confirmation")
    ns = parser.parse_args(args=args)
    autodo = ns.yes
    release_config = {}
    if ns.release is not None:
        release_config["plinder_release"] = ns.release
    if ns.release_number is not None:
        release_config["plinder_release_number"] = ns.release_number
    cfg = get_config(config={"data": release_config} if release_config else None)
    LOG.info(
        dedent(
            f"""
            Syncing {cfg.data.plinder_remote} -> {cfg.data.plinder_dir}.
            If this is the first time you are running this command, it will take a while!

            The estimated time on the progress bar may vary wildly based on varied file sizes.
            If you need to cancel this and come back to it, it will pick up where it left off.
            """
        )
    )
    for label, artifact_names, large_download in _DOWNLOAD_GROUPS:
        do_download = autodo or not large_download
        if large_download and not autodo:
            answer = input(f"Download {label}? [Y/n] ").strip().lower()
            do_download = answer in {"", "y", "yes"}
        if not do_download:
            LOG.info(f"skipping {label}; its files are fetched when requested")
            continue
        LOG.info(f"syncing {label}")
        for artifact_name in artifact_names:
            artifact_path = RELEASE_PATHS[artifact_name]
            if "{" in artifact_path:
                raise RuntimeError(
                    f"bulk download cannot resolve parameterized artifact "
                    f"{artifact_name}"
                )
            cpl.get_plinder_path(
                rel=artifact_path,
                force_progress=large_download,
            )

    t1 = time()
    total = t1 - t0
    timing = f"{total:.2f}s"
    if total > 3600:
        timing = f"{total / 3600:.2f}h"
    elif total > 60:
        timing = f"{total / 60:.2f}m"

    LOG.info(
        dedent(
            f"""
            Sync complete in {timing}!

            If you skipped large groups, plinder.core fetches their files when
            requested unless offline mode is enabled. Use
            plinder.core.download_pdb_mmcifs(...) if offline system reconstruction
            is needed.
            """
        )
    )
