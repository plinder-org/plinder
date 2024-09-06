# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
from __future__ import annotations

from functools import cached_property

import numpy as np
from PDBValidation.PDBXReader import ResidueNotFound
from PDBValidation.Residue import Residue
from PDBValidation.Validation import PDBValidation
from PDBValidation.XML import ModelledSubgroupNotFound
from pydantic import ConfigDict, Field, computed_field

from plinder.core.utils.log import setup_logger
from plinder.data.utils.annotations.utils import DocBaseModel

LOG = setup_logger(__name__)


class ResidueValidationThresholds(DocBaseModel):
    model_config = ConfigDict(ser_json_inf_nan="constants")
    min_rscc: float = 0.8
    max_rsr: float = 0.3
    min_average_occupancy: float = 1.0


def nanmean_return_nan(arr: list[float]) -> float:
    if len(arr) == 0 or all(np.isnan(arr)):
        return float(np.nan)
    return float(np.nanmean(arr))


class ResidueValidation(DocBaseModel):
    model_config = ConfigDict(ser_json_inf_nan="constants")
    altcode: str = Field(description="Alternative conformation code")
    inscode: str = Field(description="Insertion code")
    rsr: float = Field(description="Real-Space R-value")
    rsrz: float = Field(description="Real-Space R-value Z-score")
    rscc: float = Field(description="Real-Space Correlation Coefficient")
    average_occupancy: float = Field(description="Average occupancy")
    average_b_factor: float = Field(description="Average B factor")
    unknown_residue: bool = Field(description="Whether the residue is unknown")
    atom_count: int = Field(description="Number of atoms")
    unknown_atom_count: int = Field(description="Number of unknown atoms")
    heavy_atom_count: int = Field(description="Number of heavy atoms")
    num_unresolved_heavy_atoms: int = Field(
        description="Number of unresolved heavy atoms"
    )
    is_outlier: dict[str, bool] = Field(description="Whether the residue is an outlier")
    is_atom_count_consistent: bool = Field(
        description="Whether the atom count is consistent"
    )
    has_clashing_partial_occupancy_atoms: bool = Field(
        description="Whether the residue has clashing partial occupancy atoms"
    )
    alt_count: int = Field(description="Number of conformations")

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
            if (
                residue.isType("PolymerChainResidue")
                and np.abs(heavy_atom_count - heavy_atom_count_conop) <= 1
            ):
                heavy_atom_count_conop = heavy_atom_count
            num_unresolved_heavy_atoms = heavy_atom_count_conop - heavy_atom_count

            # TODO: uncomment for rerun
            #     # set -1 as value for cases where no reference was found
            #     num_unresolved_heavy_atoms = -1
            # else:
            #     num_unresolved_heavy_atoms = heavy_atom_count_conop - heavy_atom_count
            #     if (
            #         residue.isType("PolymerChainResidue")
            #         and num_unresolved_heavy_atoms >= 1
            #     ):
            #         # residues have one atom less when in polymer - adjust for it
            #         num_unresolved_heavy_atoms -= 1

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
    num_residues: int = Field(description="Number of residues in the list")
    num_processed_residues: int = Field(
        description="Number of processed residues in the list"
    )
    percent_processed_residues: float = Field(
        description="Percentage of processed residues in the list"
    )
    average_rsr: float = Field(
        description="Average Real-Space R-value across all residues in the list"
    )
    average_rsrz: float = Field(
        description="Average Real-Space R-value Z-score across all residues in the list"
    )
    average_rscc: float = Field(
        description="Average Real-Space Correlation Coefficient across all residues in the list"
    )
    average_occupancy: float = Field(
        description="Average occupancy across all residues in the list"
    )
    percent_rsr_under_threshold: float = Field(
        description="Percentage of residues with RSR under the threshold"
    )
    percent_rscc_over_threshold: float = Field(
        description="Percentage of residues with RSCC over the threshold"
    )
    percent_occupancy_over_threshold: float = Field(
        description="Percentage of residues with occupancy over the threshold"
    )
    average_b_factor: float = Field(
        description="Average B factor across all residues in the list"
    )
    unknown_residue_count: int = Field(
        description="Number of unknown residues in the list"
    )
    atom_count: int = Field(
        description="Number of atoms across all residues in the list"
    )
    heavy_atom_count: int = Field(
        description="Number of heavy atoms across all residues in the list"
    )
    num_unresolved_heavy_atoms: int = Field(
        description="Number of unresolved heavy atoms across all residues in the list"
    )
    max_alt_count: int = Field(
        description="The highest number of configurations in a single residue in the list"
    )
    percent_outliers: dict[str, float] = Field(
        description="__Percentage of outliers for each type of outlier"
    )
    # TODO: add thresholds in rerun
    # thresholds: ResidueValidationThresholds = Field(
    #     description="Thresholds used to determine if a residue is valid"
    # )

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
            average_rsr=nanmean_return_nan([residue.rsr for residue in filtered]),
            average_rsrz=nanmean_return_nan([residue.rsrz for residue in filtered]),
            average_rscc=nanmean_return_nan([residue.rscc for residue in filtered]),
            average_occupancy=nanmean_return_nan(
                [residue.average_occupancy for residue in filtered]
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
            # TODO: add thresholds in rerun
            # thresholds=residue_thresholds,
        )

    def format(self) -> dict[str, float | int | None]:
        data = {k: v for k, v in self.__dict__.items() if k != "thresholds"}
        for key in data["percent_outliers"]:
            data[f"percent_outliers_{key}"] = data["percent_outliers"][key]
        data.pop("percent_outliers")
        return {f"validation_{k}": v for k, v in data.items()}


class EntryValidation(DocBaseModel):
    model_config = ConfigDict(ser_json_inf_nan="constants")
    resolution: float = Field(description="Resolution of the PDB entry")
    rfree: float = Field(
        description="The similarity between the observed structure-factor amplitudes and those calculated from the model. Rfree should be higher than R because it is calculated using reflections not used in the refinement. See https://www.wwpdb.org/validation/XrayValidationReportHelp"
    )
    r: float = Field(
        description="The similarity between the observed structure-factor amplitudes and those calculated from the model. See https://www.wwpdb.org/validation/XrayValidationReportHelp"
    )
    clashscore: float = Field(
        description="The Molprobity Clashscore is an approximation of the overall severity of the clashes in a structure, which is defined as the number of clashes per 1000 atoms (including hydrogens). See https://www.wwpdb.org/validation/XrayValidationReportHelp"
    )
    percent_rama_outliers: float | None = Field(
        description="The percentage of Ramachandran outliers with respect to the total number of residues in the entry for which the outlier assessment is available. See https://www.wwpdb.org/validation/XrayValidationReportHelp"
    )
    percent_rota_outliers: float | None = Field(
        description="The percentage of residues with an unusual sidechain conformation with respect to the total number of residues for which the assessment is available. See https://www.wwpdb.org/validation/XrayValidationReportHelp"
    )
    data_completeness: float | None = Field(
        description="The number of expected diffraction spots is a function of data resolution and the space group. This metric describes the number of recorded reflections as a percentage of the number expected. See https://www.wwpdb.org/validation/XrayValidationReportHelp"
    )
    percent_RSRZ_outliers: float | None = Field(
        description="The percentage Real-Space R-value Z-score outliers with respect to the total number of residues for which RSRZ was computed. See https://www.wwpdb.org/validation/XrayValidationReportHelp"
    )
    atom_count: int = Field(
        description="Number of atoms in the asymmetric unit of the PDB entry"
    )
    molprobity: float | None = Field(
        description='Overall Molprobity "effective resolution", a single-score validation number based on the correlation of multiple criteria with crystallographic resolution. as described here: https://github.com/rlabduke/MolProbity/blob/6e7512e85bdea23f7ffb16e606d1f9a0abf6e5d4/cmdline/molparser.py#L662'
    )
    mean_b_factor: float = Field(
        description="Mean B-value calculated over all modelled atoms"
    )
    median_b_factor: float = Field(
        description="Median B-value calculated over all modelled atoms"
    )
    pdbx_resolution: float = Field(
        description="See https://mmcif.wwpdb.org/dictionaries/mmcif_pdbx_v50.dic/Items/_refine.ls_d_res_high.html"
    )
    pdbx_reflns_resolution: float | None = Field(
        description="See https://mmcif.wwpdb.org/dictionaries/mmcif_pdbx_v50.dic/Items/_reflns.d_resolution_high.html"
    )
    meanI_over_sigI_obs: float | None = Field(
        description="Each reflection has an intensity (I) and an uncertainty in measurement (σ(I)), thus this metric describes the signal-to-noise ratio"
    )

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

    @computed_field  # type: ignore
    @cached_property
    def r_minus_rfree(self) -> float:
        """
        The difference between r and rfree
        """
        if isinstance(self.r, float) and isinstance(self.rfree, float):
            return np.abs(self.r - self.rfree)  # type: ignore
        return np.nan

    def format(self) -> dict[str, float | None]:
        data = self.__dict__
        data["r_minus_rfree"] = self.r_minus_rfree
        return {f"validation_{k}": v for k, v in data.items()}
