# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List, Literal, Tuple, Callable

import pandas as pd

import torch
from torch.utils.data import Dataset
from torch.utils.data import DataLoader, Dataset, Sampler

from plinder.core.scores.links import query_links
from plinder.core.split.utils import get_split
from plinder.core import system

from plinder.core.structure.structure import Structure
from plinder.core.loader.transforms import StructureTransform


def structure2tensor_transform(structure: Structure) -> dict[str, torch.Tensor]:
    props: dict[str, torch.Tensor] = system.structure2tensor(
        protein_atom_coordinates=structure.protein_coords,
        protein_atom_types=structure.protein_atom_array.element,
        protein_residue_coordinates=structure.protein_coords,
        protein_residue_ids=structure.protein_atom_array.res_name,
        protein_residue_types=structure.protein_atom_array.res_id,
        ligand_mols=structure.ligand_mols,
        dtype=torch.float32,
        use_ligand_conformer=True,
    )
    return props


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
        input_structure_priority: str = "apo",
        structure_transforms: list[StructureTransform] = [],
        transform: Callable[
            [Structure], torch.Tensor | dict[str, torch.Tensor]
        ] = structure2tensor_transform,
        target_transform: Callable[
            [Structure], torch.Tensor | dict[str, torch.Tensor]
        ] = structure2tensor_transform,
        crop_equal_monomer_shapes: bool = True,
        **kwargs: Any,
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
        self._input_structure_priority = input_structure_priority
        self._structure_transforms = structure_transforms
        self._transform = transform
        self._target_transform = target_transform
        self._crop_equal_monomer_shapes = crop_equal_monomer_shapes

    def __len__(self) -> int:
        return self._num_examples

    def __getitem__(
        self, index: int
    ) -> dict[str, int | str | pd.DataFrame | dict[str, str | pd.DataFrame]]:
        if not 0 <= index < self._num_examples:
            raise IndexError(index)

        s = system.PlinderSystem(system_id=self._system_ids[index])

        if not self._file_paths_only:
            # avoid loading structure if not needed
            holo_structure, apo_structure, pred_structure = (
                s.create_masked_bound_unbound_complexes()
            )
            cropped_structures = {
                "holo": holo_structure,
                "apo": apo_structure,
                "pred": pred_structure,
            }
            input_structure = cropped_structures[self._input_structure_priority]
            target_structure = holo_structure
            for transform in self.structure_transforms:
                input_structure = transform(input_structure)
                target_structure = transform(target_structure)

            input_id = input_structure.id
            target_id = target_structure.id
            if self.transform is not None:
                input_complex = self.transform(input_structure)
            if self.target_transform is not None:
                target_complex = self.target_transform(target_structure)

        item: Dict[str, Any] = {
            "input_id": input_id,
            "target_id": target_id,
            "input_features": input_complex,
            "target": target_complex,
        }
        if self._store_file_path:
            item["path"] = s.system_cif

        return item


def get_torch_loader(
    dataset: PlinderDataset,
    batch_size: int = 2,
    shuffle: bool = True,
    sampler: "Sampler[PlinderDataset]" | None = None,
    num_workers: int = 1,
    collate_fn: Callable[[list[dict[str, Any]]], dict[str, Any]] = system.collate_batch,
    **kwargs: Any,
) -> "DataLoader[PlinderDataset]":
    return DataLoader(
        dataset,
        batch_size=batch_size,
        shuffle=shuffle,
        num_workers=num_workers,
        sampler=sampler,
        collate_fn=collate_fn,
        **kwargs,
    )


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
    cluster_column: str = "pli_qcov__70__community",
) -> pd.DataFrame:
    from torch.utils.data import WeightedRandomSampler

    if cluster_column not in split_df.columns:
        raise ValueError(f"cluster_column={cluster_column} not in split_df.columns")

    split_df = split_df[split_df.split == split]
    cluster_counts = split_df[cluster_column].value_counts().rename("cluster_count")
    split_df = split_df.merge(cluster_counts, left_on=cluster_column, right_index=True)
    split_df.reset_index(inplace=True)
    cluster_weights = 1.0 / split_df.cluster_count.values
    sampler = WeightedRandomSampler(
        weights=cluster_weights, num_samples=len(cluster_weights)
    )
    sampler_index = [i for i in sampler]
    return split_df.loc[tuple(sampler_index), ("system_id", "split")].drop_duplicates()
