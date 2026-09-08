# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
from __future__ import annotations

import biotite.structure as struc
import numpy as np
import pandas as pd
import pytest
from biotite.structure.io import pdbx
from plinder.core import PlinderSystem
from plinder.core.release import PlinderRelease
from plinder.data.annotations.cif_utils import read_mmcif_file


def test_reconstruct_linked_apo_writes_the_scored_assembly_chain(
    cif_2y4i, tmp_path, monkeypatch
):
    links_path = tmp_path / "index" / "linked_apo_structures.parquet"
    links_path.parent.mkdir()
    pd.DataFrame(
        {
            "reference_system_id": ["1abc__1__1.A__1.Z"],
            "linked_structure_id": ["2y4i_B"],
            "source_entry_id": ["2y4i"],
            "source_chain_asym_id": ["B"],
            "source_biounit_id": ["1"],
            "source_chain_instance": ["1.B"],
            "rank": [1],
        }
    ).to_parquet(links_path, index=False)

    def fetch(_release, name, **parameters):
        assert name == "linked_apo_structures"
        assert not parameters
        return links_path

    monkeypatch.setattr(PlinderRelease, "fetch", fetch)
    system = PlinderSystem(
        system_id="1abc__1__1.A__1.Z",
        reconstruction_dir=tmp_path / "reconstructed",
    )
    assert system.linked_apo_structures["linked_structure_id"].tolist() == ["2y4i_B"]

    output = system.reconstruct_linked_apo(source_mmcif=cif_2y4i)
    assert output == tmp_path / "reconstructed" / "linked_apo" / "2y4i_B.cif"
    cif_file = read_mmcif_file(output)
    block = cif_file["2y4i_B"]
    assert {
        "entry",
        "entity",
        "entity_poly",
        "entity_poly_seq",
        "chem_comp",
        "struct_asym",
        "atom_site",
    }.issubset(block)
    atoms = pdbx.get_structure(
        cif_file,
        model=1,
        use_author_fields=False,
        include_bonds=True,
    )
    assert len(np.unique(atoms.chain_id)) == 1
    assert set(block["atom_site"]["label_asym_id"].as_array(str)).issubset(
        block["struct_asym"]["id"].as_array(str)
    )
    assert set(block["atom_site"]["label_entity_id"].as_array(str)).issubset(
        block["entity"]["id"].as_array(str)
    )

    with pytest.raises(FileExistsError, match="Refusing to overwrite"):
        system.reconstruct_linked_apo(source_mmcif=cif_2y4i)

    reference = atoms.copy()
    reference.chain_id[:] = "1.A"
    reference.coord += np.array([12.0, -7.0, 3.0])
    system.__dict__["receptor_structure"] = reference
    superposed = system.superpose_linked_apo(source_mmcif=cif_2y4i)
    superposed_atoms = pdbx.get_structure(
        read_mmcif_file(superposed),
        model=1,
        use_author_fields=False,
        include_bonds=True,
    )
    reference_ca = reference[
        struc.filter_amino_acids(reference) & (reference.atom_name == "CA")
    ]
    superposed_ca = superposed_atoms[
        struc.filter_amino_acids(superposed_atoms)
        & (superposed_atoms.atom_name == "CA")
    ]
    assert struc.rmsd(reference_ca, superposed_ca) < 1e-3

    second_chain = reference.copy()
    second_chain.chain_id[:] = "1.C"
    system.__dict__["receptor_structure"] = struc.concatenate([reference, second_chain])
    with pytest.raises(ValueError, match="reference_chain is required"):
        system.superpose_linked_apo(
            source_mmcif=cif_2y4i,
            output_cif=tmp_path / "ambiguous.cif",
        )


def test_reconstruct_linked_apo_reports_when_a_system_has_no_link(
    tmp_path, monkeypatch
):
    links_path = tmp_path / "linked_apo_structures.parquet"
    pd.DataFrame(
        {
            "reference_system_id": pd.Series(dtype="string"),
            "linked_structure_id": pd.Series(dtype="string"),
            "source_entry_id": pd.Series(dtype="string"),
            "source_chain_asym_id": pd.Series(dtype="string"),
            "source_biounit_id": pd.Series(dtype="string"),
            "source_chain_instance": pd.Series(dtype="string"),
            "rank": pd.Series(dtype="int16"),
        }
    ).to_parquet(links_path, index=False)
    monkeypatch.setattr(PlinderRelease, "fetch", lambda *_args, **_kwargs: links_path)

    system = PlinderSystem(system_id="1abc__1__1.A__1.Z")
    assert system.linked_apo_structures.empty
    with pytest.raises(ValueError, match="No linked apo structure"):
        system.reconstruct_linked_apo()
