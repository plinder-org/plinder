# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np
from ost import io, mol
from ost.mol.alg.ligand_scoring_lddtpli import LDDTPLIScorer
from ost.mol.alg.ligand_scoring_scrmsd import SCRMSDScorer
from ost.mol.alg.scoring import Scorer
from posebusters import PoseBusters

from plinder.core import PlinderSystem
from plinder.core.utils.log import setup_logger

LOG = setup_logger(__name__)


@dataclass
class LigandScores:
    chain: str
    residue_number: int
    sdf_file: Path
    atom_count: int
    protein_chain_mapping: dict[str, str] | None = None
    scores: dict[str, float | None] = field(default_factory=dict)
    reference_ligand: dict[str, "LigandScores"] = field(default_factory=dict)
    posebusters: dict[str, Any] = field(default_factory=dict)


@dataclass
class ProteinScores:
    chain_mapping: dict[str, str]
    lddt: float
    bb_lddt: float
    per_chain_lddt: dict[str, float] = field(default_factory=dict)
    per_chain_bb_lddt: dict[str, float] = field(default_factory=dict)
    per_chain_coverage: dict[str, float] = field(default_factory=dict)
    score_oligo: bool = False
    oligomer_scores: dict[str, float] = field(default_factory=dict)


@dataclass
class ComplexData:
    name: str
    receptor_file: Path
    ligand_files: list[Path]
    receptor_entity: mol.EntityHandle
    ligand_views: list[mol.EntityHandle]
    num_ligands: int
    num_proteins: int

    @classmethod
    def from_plinder_system(cls, system: PlinderSystem) -> "ComplexData":
        ligand_views = []
        ligand_files = []
        for ligand in system.ligand_sdfs:
            ligand_views.append(system.ligand_views[ligand])
            ligand_files.append(Path(system.ligand_sdfs[ligand]))
        return cls(
            name=system.system_id,
            receptor_file=Path(
                system.receptor_pdb
            ),  # PoseBusters requires a PDB file, but entity is built from CIF
            ligand_files=ligand_files,
            receptor_entity=system.receptor_entity,
            ligand_views=ligand_views,
            num_ligands=system.num_ligands,
            num_proteins=system.num_proteins,
        )

    @classmethod
    def from_files(
        cls, name: str, receptor_file: Path, ligand_files: list[Path]
    ) -> "ComplexData":
        if receptor_file.suffix not in [".cif", ".pdb"]:
            raise ValueError(
                f"receptor_file must be a .cif or .pdb file, got {receptor_file}"
            )
        if receptor_file.suffix == ".cif":
            entity = io.LoadMMCIF(receptor_file.as_posix(), fault_tolerant=True)
        else:
            entity = io.LoadPDB(receptor_file.as_posix(), fault_tolerant=True)
        ligand_views = [
            io.LoadEntity(str(ligand_sdf_file), format="sdf").Select("ele != H")
            for ligand_sdf_file in ligand_files
        ]
        return cls(
            name=name,
            receptor_file=receptor_file,
            ligand_files=ligand_files,
            receptor_entity=entity,
            ligand_views=ligand_views,
            num_ligands=len(ligand_files),
            num_proteins=sum(
                1
                for x in entity.chains
                if x.type == mol.CHAINTYPE_POLY_PEPTIDE_L
                or x.type == mol.CHAINTYPE_UNKNOWN
            ),
        )


@dataclass
class ModelScores:
    reference: ComplexData
    model: ComplexData
    score_protein: bool = False
    score_posebusters: bool = False
    score_posebusters_full_report: bool = False
    posebusters_mapper: str = "bisy_rmsd"
    num_mapped_reference_ligands: int = 0
    num_mapped_model_ligands: int = 0
    ligand_scores: list[LigandScores] | None = None
    num_mapped_reference_proteins: int = 0
    num_mapped_model_proteins: int = 0
    protein_scores: ProteinScores | None = None

    @classmethod
    def from_model_files(
        cls,
        model_name: str,
        model_receptor_file: Path,
        model_ligand_files: list[Path],
        reference: PlinderSystem,
        score_protein: bool = False,
        score_posebusters: bool = False,
        score_posebusters_full_report: bool = False,
    ) -> "ModelScores":
        """
        Create a ModelScores object from a model (receptor file and a list of ligand SDF files) and a reference PlinderSystem.

        Parameters
        ----------
        model_file : Path
            The path to the model file.
        model_ligand_sdf_files : list[str | Path]
            The list of ligand SDF files.
        reference : PlinderSystem
            The reference system.
        score_protein : bool, default=False
            Whether to score the protein.
        score_posebusters : bool, default=False
            Whether to score posebusters.
        score_posebusters_full_report : bool, default=False
            Whether to score posebusters with the full report.

        Returns
        -------
        ModelScores
            The ModelScores object.
        """
        model = ComplexData.from_files(
            model_name, model_receptor_file, model_ligand_files
        )
        model_class = cls(
            reference=ComplexData.from_plinder_system(reference),
            model=model,
            score_protein=score_protein,
            score_posebusters=score_posebusters or score_posebusters_full_report,
            score_posebusters_full_report=score_posebusters_full_report,
        )
        if model_class.score_protein:
            model_class.calculate_protein_scores()
        model_class.calculate_ligand_scores()
        return model_class

    @classmethod
    def from_model_and_reference_files(
        cls,
        model_name: str,
        model_receptor_file: Path,
        model_ligand_files: list[Path],
        reference_name: str,
        reference_receptor_file: Path,
        reference_ligand_files: list[Path],
        score_protein: bool = False,
        score_posebusters: bool = False,
        score_posebusters_full_report: bool = False,
    ) -> "ModelScores":
        """
        Create a ModelScores object from both model and reference (receptor file and a list of ligand SDF files both both).

        Parameters
        ----------
        model_file : Path
            The path to the model file.
        model_ligand_sdf_files : list[str | Path]
            The list of ligand SDF files.
        reference_name : str
            The name of the reference system.
        reference_receptor_file : Path
            The path to the reference receptor file.
        reference_ligand_files : list[Path]
            The list of reference ligand SDF files.
        score_protein : bool, default=False
            Whether to score the protein.
        score_posebusters : bool, default=False
            Whether to score posebusters.
        score_posebusters_full_report : bool, default=False
            Whether to score posebusters with the full report.

        Returns
        -------
        ModelScores
            The ModelScores object.
        """
        model = ComplexData.from_files(
            model_name, model_receptor_file, model_ligand_files
        )
        reference = ComplexData.from_files(
            reference_name, reference_receptor_file, reference_ligand_files
        )
        model_class = cls(
            reference=reference,
            model=model,
            score_protein=score_protein,
            score_posebusters=score_posebusters or score_posebusters_full_report,
            score_posebusters_full_report=score_posebusters_full_report,
        )
        if model_class.score_protein:
            model_class.calculate_protein_scores()
        model_class.calculate_ligand_scores()
        return model_class

    def calculate_protein_scores(self, lddt_add_mdl_contacts: bool = False) -> None:
        scorer = Scorer(
            model=self.model.receptor_entity,
            target=self.reference.receptor_entity,
            lddt_add_mdl_contacts=lddt_add_mdl_contacts,
        )
        self.protein_scores = ProteinScores(
            chain_mapping={k: v for k, v in scorer.mapping.alns.keys()},
            lddt=scorer.lddt,
            bb_lddt=scorer.bb_lddt,
        )
        self.num_mapped_reference_proteins = len(scorer.mapping.alns)
        self.num_mapped_proteins = len(scorer.mapping.alns)
        local_lddt = scorer.local_lddt
        self.protein_scores.per_chain_lddt = {}
        bb_local_lddt = scorer.bb_local_lddt
        self.protein_scores.per_chain_bb_lddt = {}
        for reference_chain, model_chain in self.protein_scores.chain_mapping.items():
            reference_chain_length = len(
                self.reference.receptor_entity.FindChain(reference_chain).residues
            )
            self.protein_scores.per_chain_coverage[model_chain] = (
                len(local_lddt[model_chain]) / reference_chain_length
            )
            self.protein_scores.per_chain_lddt[model_chain] = (
                sum(v for v in local_lddt[model_chain].values() if v is not None)
                / reference_chain_length
            )
            self.protein_scores.per_chain_bb_lddt[model_chain] = (
                sum(v for v in bb_local_lddt[model_chain].values() if v is not None)
                / reference_chain_length
            )
        if self.reference.num_proteins > 1 and self.model.num_proteins > 1:
            self.protein_scores.score_oligo = True
            self.protein_scores.oligomer_scores = {
                "qs_global": scorer.qs_global,
                "qs_best": scorer.qs_best,
                "dockq_wave": scorer.dockq_wave,
                "dockq_ave": scorer.dockq_ave,
                # TODO: add more interface scores
            }

    def calculate_ligand_scores(self) -> None:
        pb = PoseBusters(config="dock")

        scrmsd_scorer = SCRMSDScorer(
            model=self.model.receptor_entity,
            target=self.reference.receptor_entity,
            model_ligands=self.model.ligand_views,
            target_ligands=self.reference.ligand_views,
            substructure_match=True,
            resnum_alignments=False,
        )
        lddt_pli_scorer = LDDTPLIScorer(
            model=self.model.receptor_entity,
            target=self.reference.receptor_entity,
            model_ligands=self.model.ligand_views,
            target_ligands=self.reference.ligand_views,
            substructure_match=True,
            add_mdl_contacts=True,
            resnum_alignments=False,
        )
        self.ligand_scores = []
        scrmsd_aux, lddt_pli_aux = (
            scrmsd_scorer.aux,
            lddt_pli_scorer.aux,
        )
        assigned_model = set()
        assigned_target = set()
        lddt_pli_assignment = {v: k for k, v in lddt_pli_scorer.assignment}
        scrmsd_assignment = {v: k for k, v in scrmsd_scorer.assignment}
        for m, (ligand, sdf_file) in enumerate(
            zip(self.model.ligand_views, self.model.ligand_files)
        ):
            for residue in ligand.residues:
                chain_name = residue.chain.name
                residue_number = residue.number.num
                ligand_class = LigandScores(
                    chain_name, residue_number, Path(sdf_file), len(residue.atoms)
                )
                scrmsd_aux_ligand = scrmsd_aux.get(chain_name, {}).get(
                    residue_number, {}
                )
                ligand_class.protein_chain_mapping = scrmsd_aux_ligand.get(
                    "chain_mapping", None
                )
                for name, aux, score_name, assignment in zip(
                    ["bisy_rmsd", "lddt_lp", "lddt_pli"],
                    [
                        scrmsd_aux_ligand,
                        scrmsd_aux_ligand,
                        lddt_pli_aux.get(chain_name, {}).get(residue_number, {}),
                    ],
                    ["rmsd", "lddt_lp", "lddt_pli"],
                    [scrmsd_assignment, scrmsd_assignment, lddt_pli_assignment],
                ):
                    ligand_class.scores[name] = aux.get(score_name, None)
                    if "target_ligand" in aux:
                        ligand_class.reference_ligand[name] = LigandScores(
                            aux["target_ligand"].chain.name,
                            aux["target_ligand"].number.num,
                            Path(self.reference.ligand_files[assignment[m]]),
                            len(aux["target_ligand"].atoms),
                        )
                if (
                    self.score_posebusters
                    and ligand_class.reference_ligand.get(self.posebusters_mapper, None)
                    is not None
                ):
                    result_dict = pb.bust(
                        mol_pred=ligand_class.sdf_file,
                        mol_true=ligand_class.reference_ligand[
                            self.posebusters_mapper
                        ].sdf_file,
                        mol_cond=self.reference.receptor_file,
                        full_report=self.score_posebusters_full_report,
                    ).to_dict()
                    key = (str(ligand_class.sdf_file), chain_name.split("_")[-1])
                    try:
                        ligand_class.posebusters = {
                            k: v[key] for k, v in result_dict.items()
                        }
                    except KeyError:
                        try:
                            key = (str(ligand_class.sdf_file), "mol_at_pos_0")
                            ligand_class.posebusters = {
                                k: v[key] for k, v in result_dict.items()
                            }
                        except KeyError:
                            key = (
                                str(ligand_class.sdf_file),
                                ligand_class.sdf_file.stem,
                            )
                            ligand_class.posebusters = {
                                k: v[key] for k, v in result_dict.items()
                            }
                if ligand_class.protein_chain_mapping is not None:
                    assigned_model.add(chain_name)
                    assigned_target.add(
                        ligand_class.reference_ligand["bisy_rmsd"].chain
                    )
                self.ligand_scores.append(ligand_class)
            self.num_mapped_reference_ligands = len(assigned_target)
            self.num_mapped_model_ligands = len(assigned_model)

    def get_average_ligand_scores(
        self, score_name: str
    ) -> tuple[float | None, float | None]:
        if self.ligand_scores is None:
            return None, None
        scores = []
        weights = []
        for s in self.ligand_scores:
            score = s.scores.get(score_name)
            if score is not None:
                scores.append(score)
                weights.append(s.atom_count)
        if not len(scores):
            return None, None
        average_score = np.mean(scores)
        weighted_average_score = sum(w * s for w, s in zip(weights, scores)) / sum(
            weights
        )
        return average_score, weighted_average_score

    def get_average_posebusters(self) -> dict[str, list[Any]]:
        scores: dict[str, Any] = {}
        weights: list[float] = []
        if self.ligand_scores is None:
            return scores
        for ligand_scores in self.ligand_scores:
            weights.append(ligand_scores.atom_count)
            for key, value in ligand_scores.posebusters.items():
                scores.setdefault(key, [])
                scores[key].append(value)
        avg_scores = {}
        for key, values in scores.items():
            try:
                avg_scores[f"posebusters_{key}"] = np.mean(
                    [w * v for w, v in zip(weights, values)]
                ) / np.sum(weights)
            except Exception:
                try:
                    if isinstance(values[0], bool):
                        avg_scores[f"posebusters_{key}"] = np.mean(
                            [float(val) for val in values]
                        )
                    elif isinstance(values[0], str):
                        avg_scores[f"posebusters_{key}"] = ";".join(values)
                except Exception:
                    avg_scores[f"posebusters_{key}"] = None
        return avg_scores

    def summarize_scores(self) -> dict[str, Any]:
        scores: dict[str, Any] = dict(
            model=self.model.name,
            reference=self.reference.name,
            num_reference_ligands=self.reference.num_ligands,
            num_model_ligands=self.model.num_ligands,
            num_reference_proteins=self.reference.num_proteins,
            num_model_proteins=self.model.num_proteins,
            fraction_reference_ligands_mapped=self.num_mapped_reference_ligands
            / self.reference.num_ligands,
            fraction_model_ligands_mapped=self.num_mapped_model_ligands
            / self.model.num_ligands,
        )
        score_list = ["lddt_pli", "lddt_lp", "bisy_rmsd"]
        for score_name in score_list:
            (
                scores[f"{score_name}_ave"],
                scores[f"{score_name}_wave"],
            ) = self.get_average_ligand_scores(score_name)
        if self.score_posebusters:
            scores.update(self.get_average_posebusters())
        if self.score_protein and self.protein_scores is not None:
            scores["fraction_reference_proteins_mapped"] = (
                self.num_mapped_reference_proteins / self.reference.num_proteins
            )
            scores["fraction_model_proteins_mapped"] = (
                self.num_mapped_proteins / self.model.num_proteins
            )
            scores["lddt"] = self.protein_scores.lddt
            scores["bb_lddt"] = self.protein_scores.bb_lddt
            per_chain_lddt = list(self.protein_scores.per_chain_lddt.values())
            scores["per_chain_lddt_ave"] = np.mean(per_chain_lddt)
            per_chain_bb_lddt = list(self.protein_scores.per_chain_bb_lddt.values())
            scores["per_chain_bb_lddt_ave"] = np.mean(per_chain_bb_lddt)
            if self.protein_scores.score_oligo:
                scores.update(self.protein_scores.oligomer_scores)
        return scores
