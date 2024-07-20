# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np


def fetch_pdb(pdb_id: str) -> None:
    cmd.fetch(pdb_id)
    # Takes care of edge cases where the whole protein has been modelled
    # as an alt loc.
    # select alt locs
    cmd.select("nonalts", f"{pdb_id} and poly and alt ''+'A'")
    if len(cmd.get_model("nonalts").get_residues()) == 0:
        cmd.alter(pdb_id, "alt=''")


def get_pairwise_distances(xyz1, xyz2):
    """
    Get fast pairwise Euclidean distances between coordinates.
    xyz1, xyz2 : numpy arrays of shapes (N, 3) and (M, 3)
    The result is a numpy matrix array (N, M).
    """
    A2 = np.broadcast_to([np.dot(r, r.T) for r in xyz1], (len(xyz2), len(xyz1))).T
    AB = np.dot(xyz1, xyz2.T).astype("float")
    B2 = np.broadcast_to([np.dot(r, r.T) for r in xyz2], (len(xyz1), len(xyz2)))
    return np.sqrt(A2 - 2 * AB + B2)


def get_pymol_coordinates(selection: str) -> list[list[float, float, float]]:
    # set dictionary to store variables
    pymolspace = {}
    pymolspace["xyz"] = []
    cmd.iterate_state(1, selection, "xyz.append([x, y, z])", space=pymolspace)
    return pymolspace["xyz"]


def get_cog(selection: str, ca_only: bool = True) -> tuple[float, float, float]:
    if ca_only:
        selection = f"(({selection}) and name CA)"
    xyz = get_pymol_coordinates(selection)
    return np.mean(xyz, axis=1)


def get_symmetry_contact_number(
    system_id: Path | str,
    ligand_sdf: Path | str,
    receptor_cif: Path | str,
    cutoff: float = 4,
) -> int:
    # load original pdb entry
    pdb_id = system_id.split("/")[-1].split("_")[0]
    fetch_pdb(pdb_id)

    # load ligand
    cmd.load(ligand_sdf, "ligand")

    # generate symmetry mates around the ligand
    cmd.symexp("symexp*", pdb_id, "ligand", cutoff=cutoff, segi=1)
    symmetry_mates = cmd.get_object_list("symexp*")

    # contact count
    crystal_contacts = 0

    if len(symmetry_mates):
        # detect and remove symmetry mates that overlap with receptor
        # eg. receptor incorporates more than assymetric unit to complete biounit
        cmd.load(receptor_cif, "receptor")

        # before getting coords - remove any ambiguous atoms and hydrogens
        cmd.remove("name *?")
        cmd.remove("e. h")

        cog_rec = get_cog("receptor", ca_only=True)
        xyz_lig = get_pymol_coordinates("ligand")

        for symmate in symmetry_mates:
            # get distance between CoG of CAs and compare if within tolerance
            cog_mate = get_cog(symmate, ca_only=True)
            # given some buffer for numerical error
            if np.linalg.norm(cog_mate - cog_rec) > 0.1:
                # this is a real symmetry mate that is not part of biounit
                xyz = get_pymol_coordinates(symmate)
                pair_dists = get_pairwise_distances(xyz_lig, xyz)
                print(np.shape(xyz_lig))
                print(np.shape(pair_dists))
                print(np.shape(np.min(pair_dists, axis=0)))
                # crystal_contacts as number of symmetry mate atoms within cutoff
                contact_num = np.sum(np.min(pair_dists, axis=0) < cutoff)
                crystal_contacts += contact_num
    return crystal_contacts


def main() -> None:
    # f"{PYMOL_PYTHON} -u lattice_contacts.py {system_id} {ligand_path} {receptor_path} {cutoff}"

    # use command line
    # Importing the PyMOL module will create the window.
    import pymol

    # Tell PyMOL we don't want any GUI features.
    pymol.pymol_argv = ["pymol", "-Qic"]
    # Call the function below before using any PyMOL modules.
    pymol.finish_launching()

    # import argparse
    # parser = argparse.ArgumentParser(description="Annotate and stratify a test set")
    # parser.add_argument(
    #     "--split_file",
    #     type=Path,
    #     help="Path to split file with [system_id, split] as columns",
    # )
    # parser.add_argument(
    #     "--data_dir",
    #     type=Path,
    #     help="Path to plinder data",
    # )
    # parser.add_argument(
    #     "--output_dir",
    #     type=Path,
    #     help="Path to output folder where similarity and stratification data are saved",
    # )

    # assign command line arguments
    system_id = sys.argv[1]
    ligand_sdf = sys.argv[2]
    receptor_cif = sys.argv[3]
    cutoff = int(sys.argv[4])

    # run the main code
    res = get_symmetry_contact_number(
        system_id, ligand_sdf, receptor_cif, cutoff=cutoff
    )
    print(res)


if __name__ == "__main__":
    main()
