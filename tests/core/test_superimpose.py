from zipfile import ZipFile

from biotite.sequence.io.fasta import FastaFile

from plinder.core.structure.structure import Structure


def test_superimpose_chain(read_plinder_mount, tmp_path):
    """
    Check if :func:`superimpose_chain()` can handle different scenarios.
    In all cases the superimposed structure should have the original number of atoms
    and a low RMSD to the fixed structure.
    """
    system_id_1 = "1avd__1__1.A_2.A__1.D"
    system_id_2 = "1avd__1__1.A_2.A__2.D"
    chain_id_1 = "1.A"
    with ZipFile(read_plinder_mount / "systems" / "av.zip") as archive:
        receptor_1 = archive.extract(f"{system_id_1}/receptor.cif", tmp_path)
        receptor_2 = archive.extract(f"{system_id_2}/receptor.cif", tmp_path)
        sequences_1 = archive.extract(f"{system_id_1}/sequences.fasta", tmp_path)
        sequences_2 = archive.extract(f"{system_id_2}/sequences.fasta", tmp_path)
    struct1 = Structure(
        id=system_id_1,
        protein_path=receptor_1,
        protein_sequence=dict(FastaFile.read_iter(sequences_1)),
    )
    struct2 = Structure(
        id=system_id_2,
        protein_path=receptor_2,
        protein_sequence=dict(FastaFile.read_iter(sequences_2)),
    )

    chain_1_array = struct1.protein_atom_array[
        struct1.protein_atom_array.chain_id == chain_id_1
    ]
    super_chain_1, raw_rmsd, refined_rmsd = struct1.superimpose(struct2)
    assert isinstance(super_chain_1, Structure)
    super_chain_1_array = super_chain_1.protein_atom_array[
        super_chain_1.protein_atom_array.chain_id == chain_id_1
    ]
    assert super_chain_1_array.shape == chain_1_array.shape

    # check alignment quality
    assert abs(raw_rmsd - refined_rmsd) < 0.01 or raw_rmsd > refined_rmsd
    assert refined_rmsd < 2.0
