from plinder.core.index import PlinderSystem
from plinder.core.structure.structure import Structure


def test_superimpose_chain(read_plinder_mount):
    """
    Check if :func:`superimpose_chain()` can handle different scenarios.
    In all cases the superimposed structure should have the original number of atoms
    and a low RMSD to the fixed structure.
    """
    system_id_1 = "19hc__1__1.A__1.G"
    system_id_2 = "19hc__1__1.A_1.B__1.V_1.X_1.Y"
    # system_dir_1 = read_plinder_mount / "systems" / system_id_1
    # system_dir_2 = read_plinder_mount / "systems" / system_id_2
    chain_id_1 = "1.A"
    # chain_id_2 = "1.A"
    struct1 = PlinderSystem(system_id=system_id_1).holo_structure
    struct2 = PlinderSystem(system_id=system_id_2).holo_structure

    chain_1_array = struct1.protein_atom_array[
        struct1.protein_atom_array.chain_id == chain_id_1
    ]
    # TODO: test assertions here
    # chain_2_array = sys_2.protein_atom_array[
    #     sys_2.protein_atom_array.chain_id == chain_id_2
    # ]
    super_chain_1, raw_rmsd, refined_rmsd = struct1.superimpose(struct2)
    assert isinstance(super_chain_1, Structure)
    assert super_chain_1.protein_atom_array.shape == chain_1_array.shape
    assert (raw_rmsd == 0 and refined_rmsd == 0) or raw_rmsd > refined_rmsd
    assert refined_rmsd < 2.0
