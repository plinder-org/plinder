# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
from __future__ import annotations

import json
import logging
from pathlib import Path
from typing import Any, Dict, Optional

import pandas as pd
from tqdm import tqdm

from plinder.core.utils.log import setup_logger
from plinder.data.utils.annotations.aggregate_annotations import Entry

LOG = setup_logger(__name__, log_level=logging.DEBUG)


class GetPlinderAnnotation:
    def __init__(
        self,
        mmcif_file: Path,
        validation_xml: Path,
        save_folder: Optional[Path] = None,
        neighboring_residue_threshold: float = 6.0,
        neighboring_ligand_threshold: float = 4.0,  # TODO: review @VO
        min_polymer_size: int = 10,  # TODO: review @VO,
        symmetry_mate_contact_threshold: float = 5.0,
        entry_cfg: Optional[Dict[Any, Any]] = None,
    ) -> None:
        self.mmcif_file = mmcif_file
        self.validation_xml = Path(validation_xml)
        self.save_folder = save_folder
        self.neighboring_residue_threshold = neighboring_residue_threshold
        self.neighboring_ligand_threshold = neighboring_ligand_threshold
        self.min_polymer_size = min_polymer_size
        self.symmetry_mate_contact_threshold = symmetry_mate_contact_threshold
        self.entry_cfg = entry_cfg

    def annotate(self) -> Optional[pd.DataFrame]:
        entry_cfg = dict(
            neighboring_residue_threshold=self.neighboring_residue_threshold,
            neighboring_ligand_threshold=self.neighboring_ligand_threshold,
            min_polymer_size=self.min_polymer_size,
            save_folder=self.save_folder,
            symmetry_mate_contact_threshold=self.symmetry_mate_contact_threshold,
        )
        if self.entry_cfg is not None:
            entry_cfg.update(self.entry_cfg)
        self.entry = Entry.from_cif_file(
            self.mmcif_file,
            **entry_cfg,  # type: ignore
        )
        LOG.info(f"created entry for {self.mmcif_file}")
        self.entry.set_validation(self.validation_xml, self.mmcif_file)
        if len(self.entry.systems):
            self.annotated_df = self.entry.to_df()
            if self.annotated_df.empty:
                raise ValueError(f"No ligands detected in entry {self.mmcif_file}")
            return self.annotated_df
        else:
            LOG.info(f"no entry systems for {self.mmcif_file}")
        return None


def hpc_save_batch(
    ids_file: Path,
    pdb_dir: Path,
    validation_dir: Path,
    output_file: Path,
    save_folder: Path,
) -> None:
    json_list = []
    df_list = []
    errors = []
    with open(ids_file) as f:
        for pdb_id in tqdm(f):
            pdb_id = pdb_id.strip()
            try:
                cif_file = (
                    pdb_dir
                    / pdb_id[1:3]
                    / f"pdb_0000{pdb_id}"
                    / f"pdb_0000{pdb_id}_xyz-enrich.cif.gz"
                )
                validation_file = (
                    validation_dir
                    / pdb_id[1:3]
                    / pdb_id
                    / f"{pdb_id}_validation.xml.gz"
                )
                plinder_anno = GetPlinderAnnotation(
                    cif_file,
                    validation_file,
                    save_folder=save_folder,
                )
                plinder_anno.annotate()
                json_list.append(plinder_anno.entry)
                if len(plinder_anno.entry.systems):
                    df_list.append(plinder_anno.annotated_df)
            except Exception as e:
                errors.append(f"{pdb_id}: {e}\n")
    if len(df_list):
        pd.concat(df_list).to_csv(
            output_file.with_suffix(".tsv"), index=False, sep="\t"
        )
    with open(output_file.with_suffix(".json"), "w") as f:
        f.write(json.dumps([entry.dict() for entry in json_list], indent=4))
    if len(errors):
        with open(output_file.with_suffix(".errors"), "w") as f:
            for line in errors:
                f.write(line)


def cloud_save_annotation() -> None:
    from omegaconf import OmegaConf

    from plinder.data.pipeline.config import AnnotationConfig, EntryConfig

    cfg = OmegaConf.to_container(
        OmegaConf.merge(
            {
                "mmcif_file": None,
                "validation_xml": None,
                "raise_exceptions": False,
            },
            {
                "annotation": AnnotationConfig(),
                "entry": EntryConfig(),
            },
            OmegaConf.from_cli(),
        )
    )
    assert cfg["mmcif_file"] is not None, "please pass mmcif_file=path/to/cif"
    assert cfg["validation_xml"] is not None, "please pass validation_xml=path/to/xml"
    cif = Path(cfg.pop("mmcif_file"))
    val = Path(cfg.pop("validation_xml"))
    entry_cfg = cfg.pop("entry")
    annotation_cfg = cfg.pop("annotation")
    save_folder = entry_cfg.get("save_folder")
    if save_folder is not None:
        save_folder = Path(save_folder)
        save_folder.mkdir(exist_ok=True, parents=True)
        entry_cfg["save_folder"] = save_folder
    LOG.info(f"annotating {cif}")
    gpa = GetPlinderAnnotation(
        cif,
        val,
        entry_cfg=entry_cfg,
        **annotation_cfg,
    )
    try:
        df = gpa.annotate()
        if save_folder is not None:
            pdb_id = val.stem.replace("_validation.xml", "")
            with open(f"{save_folder}/{pdb_id}.json", "w") as f:
                f.write(gpa.entry.model_dump_json(indent=4))
    except Exception as e:
        if cfg["raise_exceptions"]:
            raise e
        LOG.error(f"exception on annotating {cif}: {repr(e)}")
        if save_folder is not None:
            failure_log = Path(f"{save_folder}/failures.csv")
            pdb_id = val.stem.replace("_validation.xml", "")
            header = False
            if failure_log.is_file():
                df = pd.read_csv(failure_log, usecols=[0, 1])
                if pdb_id in df["pdb_id"]:
                    return
            else:
                header = True
            with open(f"{save_folder}/failures.csv", "a") as f:
                if header:
                    f.write("pdb_id,error\n")
                f.write(f"{pdb_id},{repr(e).replace(',', '_')}\n")


if __name__ == "__main__":
    import logging
    from sys import argv

    # TODO : this is a naive way to branch the entrypoint
    #        but will work because we in-line the entire
    #        config in the cloud_save_annotation case
    #        resulting in >> 6 CLI arguments
    if len(argv) == 6:
        logging.disable(logging.CRITICAL)
        hpc_save_batch(
            Path(argv[1]),
            Path(argv[2]),
            Path(argv[3]),
            Path(argv[4]),
            Path(argv[5]),
        )
    else:
        cloud_save_annotation()
