# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
from pathlib import Path

from plinder.core.utils.cpl import thread_pool
from plinder.core.utils.config import get_config
from plinder.core.utils.log import setup_logger

LOG = setup_logger(__name__)

def download_apo_structures(targets: list[str]) -> None:
    from plinder.data.pipeline.io import rsync_rcsb

    cfg = get_config()
    data_dir = Path(cfg.data.plinder_dir) / cfg.data.ingest

    def inner(tup: tuple[str, list[str]]) -> None: # two_char_code: str, pdb_ids: list[str]) -> None:
        two_char_code, pdb_ids = tup
        print(two_char_code, len(pdb_ids))
        if len(pdb_ids) > 20:
        # two_char_code, pdb_id = tup
            rsync_rcsb(
                kind="cif",
                two_char_code=two_char_code,
                # pdb_id=pdb_id,
                data_dir=data_dir,
            )
        else:
            for pdb_id in pdb_ids:
                rsync_rcsb(
                    kind="cif",
                    two_char_code=two_char_code,
                    pdb_id=pdb_id,
                    data_dir=data_dir,
                )

    tuples = set()
    counts = {}
    for target in targets:
        pdb_id, _ = target.split("_")
        tup = pdb_id[1:3], pdb_id
        tuples.add((pdb_id[1:3], pdb_id))
        counts.setdefault(tup[0], [])
        counts[tup[0]].append(pdb_id)
    print(counts)
    print(len(counts))
    LOG.info(f"set before dedupe: {len(tuples)}")
    tuples = sorted(tuples)
    LOG.info(f"set after dedupe: {len(tuples)}")
    LOG.info(f"downloading {len(tuples)} RCSB codes for apo structures")
    thread_pool(inner, sorted(counts.items())) # tuples)


def download_pred_structures(targets: list[str]) -> None:
    print("pred", len(targets))


def download_linked_apo_pred_structures():
    """
    Obtain the raw CIFs and PDBs for apo and pred structures
    found in the links dataset
    """

    from plinder.core.scores.links import query_links

    links = query_links(columns=["target_system", "filename"])
    groups = links.groupby("filename")
    for filename, group in groups:
        print(filename)
        if filename == "apo_links":
            download_apo_structures(group["target_system"].unique().tolist())
        elif filename == "pred_links":
            download_pred_structures(group["target_system"].unique().tolist())


if __name__ == "__main__":
    download_linked_apo_pred_structures()
