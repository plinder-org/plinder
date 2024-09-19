---
sd_hide_title: true
---

# Evaluation tools

## Evaluating docking poses across a stratified test set

The `plinder.eval` subpackage allows (1) assessing protein-ligand complex predictions against reference `plinder` systems, and
(2) correlating the performance of these predictions against the level of similarity of each test system to the corresponding training set.

The output file from running the scripts `src/plinder/eval/docking/write_scores.py` and `src/plinder/eval/docking/stratify_test_set.py` generates the same evaluation metrics as the ones we have on the public leaderboard.

The `plinder-eval` package allows

1. assessing protein-ligand complex predictions against reference `plinder` systems, and
2. correlating the performance of these predictions against the level of similarity of
each test system to the corresponding training set.

The output files from running `plinder-eval` will be used to populate the [MLSB leaderboard](https://www.mlsb.io/#challenge).

### Input

`predictions.csv` with each row representating a protein-ligand pose, and the following columns:

- `id`: An identifier for the prediction (same across different ranked poses of the same prediction)
- `reference_system_id`: `plinder` system ID to use as reference
- `receptor_file`: Path to protein CIF file. Leave blank if rigid docking, the system's receptor file will be used.
- `rank`: The rank of the pose (1-indexed)
- `confidence`: Optional score associated with the pose
- `ligand_file`: Path to the SDF file of the pose

`split.parquet` with, at a minimum, `system_id` and `split` columns mapping PLINDER systems to `train`, or `test`.

### Commands

#### Write scores

```bash
plinder_eval --prediction_file tests/test_data/eval/predictions.csv --output_dir test_eval/ --num_processes 8
```

This calculates accuracy metrics for all predicted poses compared to the reference. JSON files of each pose are stored in `test_eval/scores` and the summary file across all poses is stored in `test_eval/scores.parquet`.

The predicted pose is compared to the reference system and the following ligand scores are calculated:

- `fraction_reference_ligands_mapped`: Fraction of reference ligand chains with corresponding model chains
- `fraction_model_ligands_mapped`: Fraction of model ligand chains mapped to corresponding reference chains
- `lddt_pli_ave`: average lDDT-PLI across mapped ligands
- `lddt_pli_wave`: average lDDT-PLI across mapped ligands weighted by number of atoms
- `bisy_rmsd_ave`: average binding-site superposed symmetry-corrected RMSD across mapped ligands
- `bisy_rmsd_wave`: average binding-site superposed symmetry-corrected RMSD across mapped ligands weighted by number of atoms

If `score_protein` is set to True, the protein in `receptor_file` is compared to the system receptor file and the following scores are calculated:

- `fraction_reference_proteins_mapped`: Fraction of reference protein chains with corresponding model chains
- `fraction_model_proteins_mapped`: Fraction of model protein chains mapped to corresponding reference chains
- `lddt`: all atom lDDT
- `bb_lddt`: CA lDDT
- `per_chain_lddt_ave`: average all atom lDDT across all mapped chains
- `per_chain_lddt_wave`: average all atom lDDT across all mapped chains weighted by chain length
- `per_chain_bb_lddt_ave`: average CA lDDT across all mapped chains
- `per_chain_bb_lddt_wave`: average CA lDDT across all mapped chains weighted by chain length

For oligomeric complexes:

- `qs_global` - Global QS score
- `qs_best` - Global QS-score - only computed on aligned residues
- `dockq_ave` - Average of DockQ scores
- `dockq_wave` - Same as dockq_ave, weighted by native contacts

If `score_posebusters` is True, all posebusters checks are saved.

You can inspect the results at `test_eval/scores.parquet`

```python
>>> import pandas as pd
>>> df = pd.read_parquet("test_eval/scores.parquet")
>>> df.T
                                                   0                      1
model                              1a3b__1__1.B__1.D  1ai5__1__1.A_1.B__1.D
reference                          1a3b__1__1.B__1.D  1ai5__1__1.A_1.B__1.D
num_reference_ligands                              1                      1
num_model_ligands                                  1                      1
num_reference_proteins                             1                      2
num_model_proteins                                 1                      2
fraction_reference_ligands_mapped                1.0                    1.0
fraction_model_ligands_mapped                    1.0                    1.0
lddt_pli_ave                                 0.85815               0.510695
lddt_pli_wave                                0.85815               0.510695
bisy_rmsd_ave                               1.617184               3.665143
bisy_rmsd_wave                              1.617184               3.665143
rank                                               1                      1
```

#### Write test stratification data

(This command will not need to be run by a user, the `test_set.parquet` and `val_set.parquet` file will be provided with the split release)

```bash
plinder_stratify --split_file split.csv --output_dir test_data
```

Makes `test_data/test_set.parquet` which

- Labels the maximum similarity of each test system to the training set across all the similarity metrics
- Stratifies the test set based on training set similarity into `novel_pocket_pli`, `novel_ligand_pli`, `novel_protein`, `novel_ligand`, `novel_all` and `not_novel`
- Labels test systems with high quality.

To inspect the result of the run, do:
```python
>>> import pandas as pd
>>> df = pd.read_parquet("test_eval/test_set.parquet")
>>> df.T
                                                  0                      1
system_id                         1a3b__1__1.B__1.D  1ai5__1__1.A_1.B__1.D
pli_qcov                                        0.0                    0.0
protein_seqsim_qcov_weighted_sum                0.0                    0.0
protein_seqsim_weighted_sum                     0.0                    0.0
protein_fident_qcov_weighted_sum                0.0                    0.0
protein_fident_weighted_sum                     0.0                    0.0
protein_lddt_qcov_weighted_sum                  0.0                    0.0
protein_lddt_weighted_sum                       0.0                    0.0
protein_qcov_weighted_sum                       0.0                    0.0
pocket_fident_qcov                              0.0                    0.0
pocket_fident                                   0.0                    0.0
pocket_lddt_qcov                                0.0                    0.0
pocket_lddt                                     0.0                    0.0
pocket_qcov                                     0.0                    0.0
tanimoto_similarity_max                         0.0                    0.0
passes_quality                                False                  False
novel_pocket_pli                               True                   True
novel_ligand                                   True                   True
novel_protein                                  True                   True
novel_all                                      True                   True
not_novel                                     False                  False
>>>
```
