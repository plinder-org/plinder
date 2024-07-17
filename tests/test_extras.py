# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0

def test_get_smiles_dict(mock_alternative_datasets):
    entry_dir = mock_alternative_datasets("test")

    from plinder.data.utils.annotations.extras import get_ccd_smiles_dict

    obj = get_ccd_smiles_dict(entry_dir.parent.parent / "dbs" / "components" / "components.cif")
    assert len(obj)
