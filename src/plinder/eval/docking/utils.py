# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np
from ost import io, mol
from ost.mol.alg.ligand_scoring_lddtpli import LDDTPLIScorer
from ost.mol.alg.ligand_scoring_scrmsd import SCRMSDScorer
from ost.mol.alg.scoring import Scorer
from posebusters import PoseBusters

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
class ReferenceSystem:
    system_id: str
    receptor_cif_file: Path
    receptor_pdb_file: Path
    entity: mol.EntityHandle
    ligands: list[mol.ResidueView]
    ligand_sdf_files: dict[str, Path]
    num_ligands: int
    num_proteins: int

    @classmethod
    def from_reference_system(
        cls, system_dir: Path, reference_system: str
    ) -> "ReferenceSystem":
        cif_file = system_dir / reference_system / "receptor.cif"
        pdb_file = system_dir / reference_system / "receptor.pdb"
        entity = io.LoadMMCIF(cif_file.as_posix())
        ligand_sdf_files = {}
        ligands = []
        for chain in reference_system.split("__")[-1].split("_"):
            sdf_file = system_dir / reference_system / "ligand_files" / f"{chain}.sdf"
            ligand_sdf_files[chain] = sdf_file
            ligands.append(
                io.LoadEntity(sdf_file.as_posix(), format="sdf").Select("ele != H")
            )
        return cls(
            reference_system,
            cif_file,
            pdb_file,
            entity,
            ligands,
            ligand_sdf_files,
            len(ligands),
            len(reference_system.split("__")[-2].split("_")),
        )


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
class ModelScores:
    reference: ReferenceSystem
    system: str
    structure_file: Path
    entity: mol.EntityHandle
    ligands: list[mol.EntityHandle]
    ligand_sdf_files: list[str | Path]
    num_ligands: int
    num_proteins: int
    num_mapped_reference_ligands: int = 0
    num_mapped_ligands: int = 0
    score_posebusters: bool = False
    posebusters_mapper: str = "scrmsd"
    ligand_scores: list[LigandScores] | None = None
    rigid: bool = True
    score_protein: bool = False
    num_mapped_reference_proteins: int = 0
    num_mapped_proteins: int = 0
    protein_scores: ProteinScores | None = None

    @classmethod
    def from_files(
        cls,
        model_system: str,
        model_file: Path,
        model_ligand_sdf_files: list[str | Path],
        reference: ReferenceSystem,
        rigid: bool = True,
        score_protein: bool = False,
        score_posebusters: bool = False,
    ) -> "ModelScores":
        if model_file.suffix not in [".cif", ".pdb"]:
            raise ValueError(
                f"model_file must be a .cif or .pdb file, got {model_file}"
            )
        if model_file.suffix == ".cif":
            entity = io.LoadMMCIF(model_file.as_posix(), fault_tolerant=True)
        else:
            entity = io.LoadPDB(model_file.as_posix(), fault_tolerant=True)
        sdf_files = [
            sdf_file.as_posix() if isinstance(sdf_file, Path) else sdf_file
            for sdf_file in model_ligand_sdf_files
        ]
        ligands = [
            io.LoadEntity(ligand_sdf_file, format="sdf").Select("ele != H")
            for ligand_sdf_file in sdf_files
        ]
        model_class = cls(
            reference=reference,
            system=model_system,
            structure_file=model_file,
            entity=entity,
            ligands=ligands,
            ligand_sdf_files=model_ligand_sdf_files,
            num_ligands=len(model_ligand_sdf_files),
            num_proteins=sum(
                1 for x in entity.chains if x.type == mol.CHAINTYPE_POLY_PEPTIDE_L
            ),
            rigid=rigid,
            score_protein=score_protein,
            score_posebusters=score_posebusters,
        )
        if model_class.score_protein:
            model_class.calculate_protein_scores()
        model_class.calculate_ligand_scores()
        return model_class

    def calculate_protein_scores(self, lddt_add_mdl_contacts: bool = False) -> None:
        scorer = Scorer(
            model=self.entity,
            target=self.reference.entity,
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
                self.reference.entity.FindChain(reference_chain).residues
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
        if self.reference.num_proteins > 1 and self.num_proteins > 1:
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
            model=self.entity,
            target=self.reference.entity,
            model_ligands=self.ligands,
            target_ligands=self.reference.ligands,
            substructure_match=True,
            resnum_alignments=self.rigid,
        )
        lddt_pli_scorer = LDDTPLIScorer(
            model=self.entity,
            target=self.reference.entity,
            model_ligands=self.ligands,
            target_ligands=self.reference.ligands,
            substructure_match=True,
            add_mdl_contacts=False,
            resnum_alignments=self.rigid,
        )
        lddt_pli_scorer_amd = LDDTPLIScorer(
            model=self.entity,
            target=self.reference.entity,
            model_ligands=self.ligands,
            target_ligands=self.reference.ligands,
            substructure_match=True,
            add_mdl_contacts=True,
            resnum_alignments=self.rigid,
        )
        self.ligand_scores = []
        scrmsd_aux, lddt_pli_aux, lddt_pli_amd_aux = (
            scrmsd_scorer.aux,
            lddt_pli_scorer.aux,
            lddt_pli_scorer_amd.aux,
        )
        assigned_model = set()
        assigned_target = set()
        for ligand, sdf_file in zip(self.ligands, self.ligand_sdf_files):
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
                for name, aux, score_name in zip(
                    ["scrmsd", "lddt_lp", "lddt_pli", "lddt_pli_amd"],
                    [
                        scrmsd_aux_ligand,
                        scrmsd_aux_ligand,
                        lddt_pli_aux.get(chain_name, {}).get(residue_number, {}),
                        lddt_pli_amd_aux.get(chain_name, {}).get(residue_number, {}),
                    ],
                    ["rmsd", "lddt_lp", "lddt_pli", "lddt_pli"],
                ):
                    ligand_class.scores[name] = aux.get(score_name, None)
                    if "target_ligand" in aux:
                        ligand_class.reference_ligand[name] = LigandScores(
                            aux["target_ligand"].chain.name,
                            aux["target_ligand"].number.num,
                            self.reference.ligand_sdf_files[
                                aux["target_ligand"].chain.name.split("_")[-1]
                            ],
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
                        mol_cond=self.reference.receptor_pdb_file,
                        full_report=False,
                    ).to_dict()
                    key = (str(ligand_class.sdf_file), chain_name.split("_")[-1])
                    try:
                        ligand_class.posebusters = {
                            k: v[key] for k, v in result_dict.items()
                        }
                    except KeyError:
                        key = (str(ligand_class.sdf_file), "mol_at_pos_0")
                        ligand_class.posebusters = {
                            k: v[key] for k, v in result_dict.items()
                        }
                if ligand_class.protein_chain_mapping is not None:
                    assigned_model.add(chain_name)
                    assigned_target.add(ligand_class.reference_ligand["scrmsd"].chain)
                self.ligand_scores.append(ligand_class)
            self.num_mapped_reference_ligands = len(assigned_target)
            self.num_mapped_ligands = len(assigned_model)

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

    def get_passing_posebusters(self) -> dict[str, bool]:
        scores = defaultdict(list)
        if self.ligand_scores is None:
            return {}
        for ligand_scores in self.ligand_scores:
            for k in ligand_scores.posebusters:
                scores[k].append(ligand_scores.posebusters[k])
        return {f"posebusters_{k}": all(v) for k, v in scores.items()}

    def summarize_scores(self) -> dict[str, Any]:
        scores = dict(
            model=self.system,
            reference=self.reference.system_id,
            num_reference_ligands=self.reference.num_ligands,
            num_model_ligands=self.num_ligands,
            num_reference_proteins=self.reference.num_proteins,
            num_model_proteins=self.num_proteins,
            fraction_reference_ligands_mapped=self.num_mapped_reference_ligands
            / self.reference.num_ligands,
            fraction_model_ligands_mapped=self.num_mapped_ligands / self.num_ligands,
        )
        score_list = ["lddt_pli", "lddt_pli_amd", "scrmsd"]
        if self.score_protein:
            score_list.append("lddt_lp")
        for score_name in score_list:
            (
                scores[f"{score_name}_ave"],
                scores[f"{score_name}_wave"],
            ) = self.get_average_ligand_scores(score_name)
        if self.score_posebusters:
            scores.update(self.get_passing_posebusters())
        if self.score_protein and self.protein_scores is not None:
            scores["fraction_reference_proteins_mapped"] = (
                self.num_mapped_reference_proteins / self.reference.num_proteins
            )
            scores["fraction_model_proteins_mapped"] = (
                self.num_mapped_proteins / self.num_proteins
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
