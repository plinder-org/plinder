# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
from __future__ import annotations

import numpy as np
from PDBValidation.PDBXReader import ResidueNotFound
from PDBValidation.Residue import Residue
from PDBValidation.Validation import PDBValidation
from PDBValidation.XML import ModelledSubgroupNotFound
from pydantic import ConfigDict

from plinder.core.utils.log import setup_logger
from plinder.data.structure.models import DocBaseModel

LOG = setup_logger(__name__)


class ResidueValidationThresholds(DocBaseModel):
    model_config = ConfigDict(ser_json_inf_nan="constants")
    min_rscc: float = 0.8
    max_rsr: float = 0.3
    min_average_occupancy: float = 1.0


class ResidueValidation(DocBaseModel):
    model_config = ConfigDict(ser_json_inf_nan="constants")
    altcode: str
    inscode: str
    rsr: float
    rsrz: float
    rscc: float
    average_occupancy: float
    average_b_factor: float
    unknown_residue: bool
    atom_count: int
    unknown_atom_count: int
    heavy_atom_count: int
    num_unresolved_heavy_atoms: int
    is_outlier: dict[str, bool]
    is_atom_count_consistent: bool
    has_clashing_partial_occupancy_atoms: bool
    alt_count: int

    @classmethod
    def from_residue(
        cls, chain: str, resnum: str | int, entity: str, doc: PDBValidation
    ) -> ResidueValidation | None:
        try:
            residue_with_alts = Residue.CreateFromMmCIFPosition(
                chain, resnum, str(entity), doc
            )
        except (ResidueNotFound, ModelledSubgroupNotFound) as e:
            LOG.warning(
                f"from_residue: Error creating residue {chain} {resnum} {entity}: {e}"
            )
            return None
        alts = residue_with_alts.listAlt()
        if "." in alts:
            alt = "."
        else:
            alt = list(alts)[0]

        try:
            residue = residue_with_alts.getAlt(alt)
            heavy_atom_count = residue.countAtomsHeavyPDBX()
            heavy_atom_count_conop = residue.countAtomsHeavyConop()
            if heavy_atom_count_conop is None:
                heavy_atom_count_conop = 0
                # set -1 as value for cases where no reference was found
                num_unresolved_heavy_atoms = -1
            else:
                num_unresolved_heavy_atoms = heavy_atom_count_conop - heavy_atom_count
                if (
                    residue.isType("PolymerChainResidue")
                    and num_unresolved_heavy_atoms >= 1
                ):
                    # residues have one atom less when in polymer - adjust for it
                    num_unresolved_heavy_atoms -= 1
            return cls(
                altcode=residue.getAltCode(),
                inscode=residue.getInsertionCode(),
                alt_count=len(alts),
                rsr=residue.getRSR(),
                rsrz=residue.getRSRZ(),
                rscc=residue.getRSCC(),
                average_occupancy=np.nanmean(residue.getOccupancies()),
                average_b_factor=np.nanmean(residue.getBFactors()),
                unknown_residue=residue.getResName() in PDBValidation.unknown_resnames,
                atom_count=residue.countAtomsPDBX(),
                heavy_atom_count=heavy_atom_count,
                num_unresolved_heavy_atoms=num_unresolved_heavy_atoms,
                unknown_atom_count=residue.countUnknownAtoms(),
                is_outlier={
                    s: residue.isOutlier(s)
                    for s in ["geometry", "density", "chirality", "clashes"]
                },
                is_atom_count_consistent=residue.isAtomCountConsistent(),
                has_clashing_partial_occupancy_atoms=residue.hasClashingPartialOccupancyAtoms(),
            )
        except Exception as e:
            LOG.error(
                f"from_residue: Error getting residue {chain} {resnum} {entity}: {e}"
            )
            return None


class ResidueListValidation(DocBaseModel):
    model_config = ConfigDict(ser_json_inf_nan="constants")
    num_residues: int
    num_processed_residues: int
    percent_processed_residues: float
    average_rsr: float
    average_rsrz: float
    average_rscc: float
    average_occupancy: float
    percent_rsr_under_threshold: float
    percent_rscc_over_threshold: float
    percent_occupancy_over_threshold: float
    average_b_factor: float
    unknown_residue_count: int
    atom_count: int
    heavy_atom_count: int
    num_unresolved_heavy_atoms: int
    max_alt_count: int
    percent_outliers: dict[str, float]

    @classmethod
    def from_residues(
        cls,
        residues: list[ResidueValidation],
        residue_thresholds: ResidueValidationThresholds,
    ) -> ResidueListValidation | None:
        total_residues = len(residues)
        filtered = [residue for residue in residues if residue is not None]
        num_processed = len(filtered)
        if num_processed == 0:
            return None
        return cls(
            num_residues=total_residues,
            num_processed_residues=num_processed,
            percent_processed_residues=100 * num_processed / total_residues,
            average_rsr=float(np.nanmean([residue.rsr for residue in filtered])),
            average_rsrz=float(np.nanmean([residue.rsrz for residue in filtered])),
            average_rscc=float(np.nanmean([residue.rscc for residue in filtered])),
            average_occupancy=float(
                np.nanmean([residue.average_occupancy for residue in filtered])
            ),
            percent_rsr_under_threshold=100
            * sum([residue.rsr <= residue_thresholds.max_rsr for residue in filtered])
            / total_residues,
            percent_rscc_over_threshold=100
            * sum([residue.rscc > residue_thresholds.min_rscc for residue in filtered])
            / total_residues,
            percent_occupancy_over_threshold=100
            * sum(
                [
                    residue.average_occupancy
                    >= residue_thresholds.min_average_occupancy
                    for residue in filtered
                ]
            )
            / total_residues,
            average_b_factor=float(
                np.mean([residue.average_b_factor for residue in filtered])
            ),
            unknown_residue_count=sum(
                [residue.unknown_residue for residue in filtered]
            ),
            atom_count=sum([residue.atom_count for residue in filtered]),
            heavy_atom_count=sum([residue.heavy_atom_count for residue in filtered]),
            num_unresolved_heavy_atoms=sum(
                [residue.num_unresolved_heavy_atoms for residue in filtered]
            ),
            max_alt_count=max([residue.alt_count for residue in filtered]),
            percent_outliers={
                s: 100
                * sum([residue.is_outlier[s] for residue in filtered])
                / total_residues
                for s in ["geometry", "density", "chirality", "clashes"]
            },
        )

    def to_dict(self) -> dict[str, float]:
        # TODO : compare with self.model_dump
        data = {k: v for k, v in self.__dict__.items()}
        for key in data["percent_outliers"]:
            data[f"percent_outliers_{key}"] = data["percent_outliers"][key]
        data.pop("percent_outliers")
        return {f"validation_{k}": v for k, v in data.items()}


class EntryValidation(DocBaseModel):
    model_config = ConfigDict(ser_json_inf_nan="constants")
    resolution: float
    rfree: float
    r: float
    clashscore: float
    percent_rama_outliers: float | None
    percent_rota_outliers: float | None
    data_completeness: float | None
    percent_RSRZ_outliers: float | None
    atom_count: int
    molprobity: float | None
    mean_b_factor: float
    median_b_factor: float
    pdbx_resolution: float
    pdbx_reflns_resolution: float | None
    meanI_over_sigI_obs: float | None

    @classmethod
    def from_entry(cls, doc: PDBValidation) -> EntryValidation:
        xml = doc.getValidationXML()
        entry = xml.getEntry()
        rfree = entry.get("PDB-Rfree")
        if rfree == "NotAvailable":
            rfree = np.nan
        try:
            reflns = doc.getReflectionsResolution()
        except KeyError:
            reflns = None
        meanI_over_sigI_obs = doc.getMeanIOverSigIObs()
        if meanI_over_sigI_obs == "?":
            meanI_over_sigI_obs = None
        return cls(
            resolution=entry.get("PDB-resolution"),
            rfree=rfree,
            r=entry.get("PDB-R"),
            clashscore=entry.get("clashscore"),
            percent_rama_outliers=entry.get("percent-rama-outliers"),
            percent_rota_outliers=entry.get("percent-rota-outliers"),
            data_completeness=entry.get("DataCompleteness"),
            percent_RSRZ_outliers=entry.get("percent-RSRZ-outliers"),
            atom_count=doc.countAtoms(),
            molprobity=doc.calcMolProbityOverallScore(),
            mean_b_factor=doc.calcMeanStructureBFactor(),
            median_b_factor=doc.calcMedianStructureBFactor(),
            pdbx_resolution=doc.getResolution(),
            pdbx_reflns_resolution=reflns,
            meanI_over_sigI_obs=meanI_over_sigI_obs,
        )

    @property
    def r_minus_rfree(self) -> float:
        if isinstance(self.r, float) and isinstance(self.rfree, float):
            return np.abs(self.r - self.rfree)  # type: ignore
        return np.nan

    def to_dict(self) -> dict[str, float | None]:
        data = self.__dict__
        data["r_minus_rfree"] = self.r_minus_rfree
        return {f"validation_{k}": v for k, v in data.items()}
