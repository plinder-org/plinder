# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
import pandas as pd
from plinder.data.final_structure_qc import run_all_checks


def test_final_structure_checks(
    mini_system_dir, target_structure_validation_file, mini_all_json
):
    df = run_all_checks(mini_system_dir, mini_all_json)
    df["ligand_molvs_validation"] = df.ligand_molvs_validation.astype("str")
    df["ligand_rdkit_validation"] = df.ligand_rdkit_validation.astype("str")
    target_df = pd.read_csv(target_structure_validation_file, sep="\t")
    pd.testing.assert_frame_equal(df, target_df.fillna(""))
