# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
from __future__ import annotations

import sys
from pathlib import Path
from typing import Any, Dict, List, Tuple

import numpy as np
from numpy.typing import NDArray


def fetch_pdb(pdb_id: str) -> None:
    # TODO: cmd.load(mmcif_file) - get from nextgen source
    cmd.fetch(pdb_id)
    # Takes care of edge cases where the whole protein has been modelled
    # as an alt loc.
    # select alt locs
    cmd.select("nonalts", f"{pdb_id} and poly and alt ''+'A'")
    if len(cmd.get_model("nonalts").get_residues()) == 0:
        cmd.alter(pdb_id, "alt=''")


def get_pairwise_distances(
    xyz1: NDArray[np.float64], xyz2: NDArray[np.float64]
) -> NDArray[np.float64]:
    """
    Get fast pairwise Euclidean distances between coordinates.
    xyz1, xyz2 : numpy arrays of shapes (N, 3) and (M, 3)
    The result is a numpy matrix array (N, M).
    """
    A2 = np.broadcast_to([np.dot(r, r.T) for r in xyz1], (len(xyz2), len(xyz1))).T
    AB = np.dot(xyz1, xyz2.T).astype("float")
    B2 = np.broadcast_to([np.dot(r, r.T) for r in xyz2], (len(xyz1), len(xyz2)))
    return np.sqrt(A2 - 2 * AB + B2)


def get_pymol_coordinates(selection: str) -> NDArray[Any]:
    # set dictionary to store variables
    pymolspace: Dict[str, List[Any]] = {}
    pymolspace["xyz"] = []
    cmd.iterate_state(1, selection, "xyz.append([x, y, z])", space=pymolspace)
    return np.array(pymolspace["xyz"])


def get_cog(selection: str, ca_only: bool = True) -> NDArray[np.float64]:
    if ca_only:
        selection = f"(({selection}) and name CA)"
    xyz = get_pymol_coordinates(selection)
    return np.mean(xyz, axis=0)


def get_symmetry_contact_number(
    system_id: str,
    ligand_sdf: Path | str,
    receptor_cif: Path | str,
    cutoff: float = 4,
) -> Tuple[int, int]:
    # load original pdb entry
    pdb_id = system_id.split("_")[0]
    fetch_pdb(pdb_id)

    # load ligand
    cmd.load(ligand_sdf, "ligand")

    # generate symmetry mates around the ligand
    cmd.symexp("symexp*", pdb_id, "ligand", cutoff=cutoff, segi=1)
    symmetry_mates = cmd.get_object_list("symexp*")
    # TODO: test for cases when symmetry mates are not & polymer

    # contact count
    num_symmetry_mates_contacts = 0
    num_lig_contacts = 0

    if len(symmetry_mates):
        # detect and remove symmetry mates that overlap with receptor
        # eg. receptor incorporates more than assymetric unit to complete biounit
        cmd.load(receptor_cif, "receptor")
        cmd.split_chains("receptor")
        receptor_chains = cmd.get_object_list("receptor_*")

        # before getting coords - remove any ambiguous atoms and hydrogens
        cmd.remove("name *?")
        cmd.remove("e. h")

        receptor_cogs = []
        for rec_chain in receptor_chains:
            receptor_cogs.append(get_cog(f"({rec_chain} & polymer)", ca_only=True))
        xyz_lig = get_pymol_coordinates("ligand")
        # store ligand atom contacts as array in case of multiple symmetry mates in contact
        lig_atoms_in_crystal_contact = np.zeros(len(xyz_lig))

        for symmate in symmetry_mates:
            # get distance between CoG of CAs and compare if within tolerance
            cog_mate = get_cog(f"({symmate} & polymer)", ca_only=True)
            # given some buffer for numerical error
            # only consider symmetry mates that are not part of receptor
            if sum([np.linalg.norm(cog_mate - cog_rec) < 0.5 for cog_rec in receptor_cogs]) > 0:
                # this is part of biounit!
                cmd.delete(symmate)
                continue
            # this is a real symmetry mate that is not part of biounit
            xyz = get_pymol_coordinates(symmate)
            pair_contacts = get_pairwise_distances(xyz_lig, xyz) < cutoff
            lig_atoms_in_crystal_contact += np.max(pair_contacts, axis=1)
            # print(lig_atoms_in_crystal_contact)
            # print(np.shape(pair_contacts))
            # crystal_contacts as number of symmetry mate atoms within cutoff
            num_symmetry_mates_contacts += np.sum(np.max(pair_contacts, axis=0))
        # count each ligand atom only once!
        num_lig_contacts = np.sum(lig_atoms_in_crystal_contact > 0)
        # TODO: test for multiple symmetry mates in contact
        # print(num_lig_contacts, num_symmetry_mates_contacts)
    return num_lig_contacts, num_symmetry_mates_contacts


def main() -> None:
    # assign command line arguments
    system_id = sys.argv[1]
    ligand_sdf = sys.argv[2]
    receptor_cif = sys.argv[3]
    output_file = sys.argv[4]
    cutoff = float(sys.argv[5])

    # run the main code
    num_cont_lig, num_cont_symm = get_symmetry_contact_number(
        system_id, ligand_sdf, receptor_cif, cutoff=cutoff
    )
    # print(num_cont_lig, num_cont_symm)
    # print(','.join([str(num_cont_lig), str(num_cont_symm)]))
    with open(output_file, "w") as out:
        out.write(",".join([str(num_cont_lig), str(num_cont_symm)]))
    # save file to check visually
    cmd.select('lattice_contact', 'symexp* within 5 of ligand')
    cmd.show("spheres", "lattice_contact")
    cmd.save(output_file + '.pse')


if __name__ == "__main__":
    # use command line
    # f"{PYMOL_PYTHON} -u pymol_lattice_contacts_script.py {system_id} {ligand_path} {receptor_path} {output_file} {cutoff}"
    # python pymol_lattice_contacts_script.py 1az2__1__1.A__1.B_1.C ~/plinder_local_data/v2/az/1az2__1__1.A__1.B_1.C/ligand_files/1.B.sdf ~/plinder_local_data/v2/az/1az2__1__1.A__1.B_1.C/receptor.cif out.csv 6

    import pymol
    import pymol.cmd as cmd

    # Importing the PyMOL module will create the window.
    # Tell PyMOL we don't want any GUI features.
    pymol.pymol_argv = ["pymol", "-Qic"]
    # Call the function below before using any PyMOL modules.
    pymol.finish_launching()
    main()
