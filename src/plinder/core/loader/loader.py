# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List, Literal, Tuple

import atom3d.util.formats as fo
import pandas as pd
from torch.utils.data import Dataset

from plinder.core.scores.links import query_links
from plinder.core.split.utils import get_split
from plinder.core.system import system


class PlinderDataset(Dataset):  # type: ignore
    """
    Creates a dataset from plinder systems

    Parameters
    ----------
    df : pd.DataFrame | None
        the split to use
    split : str
        the split to sample from
    file_with_system_ids : str | Path
        path to a file containing a list of system ids (default: full index)
    store_file_path : bool, default=True
        if True, include the file path of the source structures in the dataset
    num_alternative_structures : int, default=0
        if available, load up to this number of alternative structures (apo and pred)
    """

    def __init__(
        self,
        df: pd.DataFrame | None = None,
        split: str = "train",
        split_parquet_path: str | Path | None = None,
        store_file_path: bool = True,
        num_alternative_structures: int = 0,
        file_paths_only: bool = False,
    ):
        if df is None:
            if split_parquet_path is None:
                df = get_split()
            else:
                df = pd.read_parquet(split_parquet_path)
        self._system_ids = df.loc[df["split"] == split, "system_id"].to_list()
        self._num_examples = len(self._system_ids)
        self._store_file_path = store_file_path or file_paths_only
        self._links = None
        if num_alternative_structures > 0:
            self._links = query_links().groupby("reference_system_id")
        self.num_alternative_structures = num_alternative_structures
        self._file_paths_only = file_paths_only

    def __len__(self) -> int:
        return self._num_examples

    def __getitem__(
        self, index: int
    ) -> dict[str, int | str | pd.DataFrame | dict[str, str | pd.DataFrame]]:
        if not 0 <= index < self._num_examples:
            raise IndexError(index)

        s = system.PlinderSystem(system_id=self._system_ids[index])

        if self._file_paths_only:
            # avoid loading structure if not needed
            structure_df = None
        else:
            structure_df = fo.bp_to_df(fo.read_any(s.system_cif))

        item: Dict[str, Any] = {
            "id": index,
            "system_id": s.system_id,
            "df": structure_df,
            "alternative_structures": {},
        }
        if self._store_file_path:
            item["path"] = s.system_cif

        if self._links is not None:
            alternatives: dict[str, pd.DataFrame | str | None] = {}
            try:
                links = self._links.get_group(s.system_id)
                if not links.empty:
                    alts = links.groupby("kind").head(self.num_alternative_structures)
                    # TODO: make better stagger between the kinds!
                    count = 0
                    for kind, link_id in zip(alts["kind"], alts["id"]):
                        structure = s.get_linked_structure(
                            link_kind=str(kind), link_id=link_id
                        )
                        df_key = f"{kind}_{link_id}_df"
                        if self._file_paths_only:
                            alternatives[df_key] = None
                        else:
                            alternatives[df_key] = fo.bp_to_df(fo.read_any(structure))
                        if self._store_file_path:
                            alternatives[f"{kind}_{link_id}_path"] = structure
                        count += 1
                        # if we have enough alternatives, return them
                        if count >= self.num_alternative_structures:
                            item["alternative_structures"] = alternatives
                            return item
            except KeyError:
                pass
            item["alternative_structures"] = alternatives
        return item


def get_model_input_files(
    split_df: pd.DataFrame,
    split: Literal["train", "val", "test"] = "train",
    max_num_sample: int = 10,
    num_alternative_structures: int = 1,
) -> List[Tuple[Path, str, List[Any]]]:
    model_inputs = []

    smiles_in_df = True
    if "ligand_rdkit_canonical_smiles" not in split_df.columns:
        smiles_in_df = False
        from rdkit import Chem

    dataset = PlinderDataset(
        df=split_df,
        split=split,
        num_alternative_structures=num_alternative_structures,
        file_paths_only=True,
    )

    count = 0
    for data in dataset:
        system_dir = Path(data["path"]).parent
        protein_fasta_filepath = system_dir / "sequences.fasta"

        if num_alternative_structures:
            alt_structure_filepaths = [
                val
                for key, val in data["alternative_structures"].items()
                if key.endswith("_path")
            ]
        else:
            alt_structure_filepaths = []

        if smiles_in_df:
            smiles_array = split_df.loc[
                split_df["system_id"] == data["system_id"],
                "ligand_rdkit_canonical_smiles",
            ].values
        else:
            smiles_array = [
                Chem.MolToSmiles(next(Chem.SDMolSupplier(ligand_sdf)))
                for ligand_sdf in (system_dir / "ligand_files/").glob("*sdf")
            ]
        # in case there's more than one!
        smiles = ".".join(smiles_array)
        model_inputs.append((protein_fasta_filepath, smiles, alt_structure_filepaths))
        count += 1
        if count >= max_num_sample:
            break
    return model_inputs


# Example sampler function, this is for demonstration purposes,
# users are advised to write a sampler that suit their need
def get_diversity_samples(
    split_df: pd.DataFrame,
    split: Literal["train", "val", "test"] = "train",
    cluster_tag: None | str = None,
    metric: Literal[
        "pli_qcov",
        "pli_unique_qcov",
        "pocket_fident",
        "pocket_fident_qcov",
        "pocket_lddt",
        "pocket_lddt_qcov",
        "pocket_qcov",
        "protein_fident_max",
        "protein_fident_qcov_max",
        "protein_fident_qcov_weighted_max",
        "protein_fident_qcov_weighted_sum",
        "protein_fident_weighted_max",
        "protein_fident_weighted_sum",
        "protein_lddt_max",
        "protein_lddt_qcov_max",
        "protein_lddt_qcov_weighted_max",
        "protein_lddt_qcov_weighted_sum",
        "protein_lddt_weighted_max",
        "protein_lddt_weighted_sum",
        "protein_qcov_weighted_sum",
        "protein_seqsim_max",
        "protein_seqsim_qcov_max",
        "protein_seqsim_qcov_weighted_max",
        "protein_seqsim_qcov_weighted_sum",
        "protein_seqsim_weighted_max",
        "protein_seqsim_weighted_sum",
    ] = "pli_qcov",
    directed: bool = True,
    threshold: int = 100,
    cluster_type: Literal["component", "communities"] = "component",
) -> pd.DataFrame:
    from torch.utils.data import WeightedRandomSampler

    if not cluster_tag:
        direction_tag = "strong" if directed else "weak"
        cluster_tag = f"{metric}__{threshold}__{direction_tag}__{cluster_type}"

    split_df = split_df[split_df.split == split]
    cluster_counts = split_df[cluster_tag].value_counts().rename("cluster_count")
    split_df = split_df.merge(cluster_counts, left_on=cluster_tag, right_index=True)
    split_df.reset_index(inplace=True)
    cluster_weights = 1.0 / split_df.cluster_count.values
    sampler = WeightedRandomSampler(
        weights=cluster_weights, num_samples=len(cluster_weights)
    )
    sampler_index = [i for i in sampler]
    return split_df.loc[tuple(sampler_index), ("system_id", "split")].drop_duplicates()
