# Test data sources

## boltz_8c3u_input_model_0.cif + boltz_8c3u_input.yaml

Boltz prediction for PDB 8c3u (protein-ligand complex). The CIF is the
predicted structure output; the YAML is the exact input config given to
Boltz (protein sequence + ligand SMILES). Tests parse the SMILES from the
YAML so there's no duplicate source of truth.

- CIF source: https://github.com/plinder-org/runs-n-poses/blob/main/examples/outputs/boltz/8c3u__1__1.A__1.C/1372115236/boltz_results_input/predictions/input/input_model_0.cif
- YAML source: https://github.com/plinder-org/runs-n-poses/blob/main/examples/inputs/boltz/8c3u__1__1.A__1.C/input.yaml
- License: Apache-2.0 (plinder-org/runs-n-poses repository)
- Used to test custom CIF processing when `_chem_comp_bond` is absent.
