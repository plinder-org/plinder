# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
import numpy as np
from plinder.data.smallmolecules import (
    NumHAcceptors,
    NumHDonors,
    containC,
    fcsp3,
    get_ecfp_fingerprint,
    mw,
    n_ring,
    n_ro_bonds,
    neutralize_mol,
    nha,
    rdkit_valid,
    smil2inchikey,
    smil2nochiral_nocharge,
)


def test_smallmolecules():
    smiles = "CC(=O)OC1=CC=CC=C1C(=O)O"
    assert nha(smiles) == 13
    assert mw(smiles) == 180.04
    assert n_ring(smiles) == 1
    assert n_ro_bonds(smiles) == 2
    assert fcsp3(smiles) == 0.11
    assert NumHAcceptors(smiles) == 3
    assert NumHDonors(smiles) == 1
    assert smil2inchikey(smiles) == "BSYNRYMUTXBXSQ-UHFFFAOYSA-N"
    assert neutralize_mol(smiles) == None
    assert smil2nochiral_nocharge(smiles) == "CC(=O)Oc1ccccc1C(=O)O"
    assert rdkit_valid(smiles) == True
    assert containC(smiles) == True
    assert np.allclose(
        np.array(
            [
                389,
                456,
                541,
                589,
                650,
                695,
                740,
                807,
                909,
                1017,
                1027,
                1035,
                1047,
                1057,
                1088,
                1106,
                1199,
                1380,
                1410,
                1447,
                1467,
                1468,
                1616,
                1729,
                1750,
                1775,
                1865,
                1873,
                1917,
                1970,
                1991,
            ]
        ),
        np.nonzero(get_ecfp_fingerprint(smiles)),
    )
