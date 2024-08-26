from pathlib import Path

import pytest

from plinder.core import system

@pytest.mark.parametrize("system_id", [
    "19hc__1__1.A_1.B__1.D_1.L_1.Q_1.S_1.U",
    "19hc__1__1.A_1.B__1.E_1.F_1.H_1.J_1.O",
    "19hc__1__1.A_1.B__1.K_1.M_1.N",
    "19hc__1__1.A_1.B__1.V_1.X_1.Y",
    "19hc__1__1.A__1.G",
    "19hc__1__1.A__1.I",
])
def test_plinder_system(system_id, read_plinder_mount):
    system.PlinderSystem(system_id=system_id)


@pytest.mark.parametrize("system_id", [
    "19hc__1__1.A__1.C",
    "19hc__1__1.B__1.P",
])
def test_plinder_system_fails(system_id, read_plinder_mount):
    with pytest.raises(ValueError):
        system.PlinderSystem(system_id=system_id).system


def test_plinder_system_system_files(read_plinder_mount):
    system_id = "19hc__1__1.A_1.B__1.V_1.X_1.Y"
    s = system.PlinderSystem(system_id=system_id)
    assert len(s.structures) == 10
    assert len(s.ligands) == 3
    assert len(s.system_cif)
    assert len(s.receptor_cif)
    assert len(s.receptor_pdb)
    assert len(s.sequences)
    assert len(s.chain_mapping)
    assert len(s.water_mapping)
    assert Path(s.system_cif).is_file()
    assert Path(s.receptor_cif).is_file()
    assert Path(s.receptor_pdb).is_file()
    assert Path(s.sequences).is_file()
    assert isinstance(s.chain_mapping, dict)
    assert isinstance(s.water_mapping, dict)
