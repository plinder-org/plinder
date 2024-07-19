# Evaluating docking poses across a stratified test set

The `plinder-eval` package allows (1) assessing protein-ligand complex predictions against reference `plinder` systems, and
(2) correlating the performance of these predictions against the level of similarity of each test system to the corresponding training set.

The output file from running `plinder-eval` can be used as a submission to the leaderboard (coming soon).

## Input
`predictions.csv` with each row representating a protein-ligand pose, and the following columns:
- `id`: An identifier for the prediction (same across different ranked poses of the same prediction)
- `reference_system_id`: `plinder` system ID to use as reference
- `receptor_file`: Path to protein CIF file. Leave blank if rigid docking, the system's receptor file will be used.
- `rank`: The rank of the pose (1-indexed)
- `confidence`: Optional score associated with the pose
- `ligand_file`: Path to the SDF file of the pose

`split.csv` with `system_id` and `split` columns mapping PLINDER systems to `train`, or `test`.

## Commands

### Write scores
```bash
python src/plinder/eval/docking/write_scores.py --prediction_file predictions.csv --data_dir PLINDER_DATA_DIR --output_dir scores --num_processes 64
```
This calculates accuracy metrics for all predicted poses compared to the reference. JSON files of each pose are stored in `scores/scores` and the summary file across all poses is stored in `scores.parquet`.

The predicted pose is compared to the reference system and the following ligand scores are calculated:
- `fraction_reference_ligands_mapped`: Fraction of reference ligand chains with corresponding model chains
- `fraction_model_ligands_mapped`: Fraction of model ligand chains mapped to corresponding reference chains
- `lddt_pli_ave`: average lDDT-PLI across mapped ligands
- `lddt_pli_wave`: average lDDT-PLI across mapped ligands weighted by number of atoms
- `lddt_pli_amd_ave`: average lDDT-PLI including penalising added model contacts across mapped ligands
- `lddt_pli_wave`: average lDDT-PLI including penalising added model contacts across mapped ligands weighted by number of atoms
- `scrmsd_ave`: average symmetry-corrected RMSD across mapped ligands
- `scrmsd_wave`: average symmetry-corrected RMSD across mapped ligands weighted by number of atoms

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


### Write test stratification data
(This command will not need to be run by a user, the `test_set.parquet` and `val_set.parquet` file will be provided with the split release)

```bash
python src/plinder/eval/docking/stratify_test_set.py --split_file split.csv --data_dir PLINDER_DATA_DIR --output_dir test_data --num_processes 16
```

Makes `test_data/test_set.parquet` which
- Labels the maximum similarity of each test system to the training set across all the similarity metrics
- Stratifies the test set based on training set similarity into `novel_pocket_pli`, `novel_pocket_ligand`, `novel_protein`, `novel_all`, and `not_novel`
- Labels test systems with high quality


### Write evaluation results

```bash
python src/plinder/eval/docking/make_plots.py --score_file scores/scores.parquet --data_file test_data/test_set.parquet --output_dir results
```

Writes out results.csv and plots of performance as a function of training set similarity across different similarity metrics.
