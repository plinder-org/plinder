![plinder](https://github.com/user-attachments/assets/05088c51-36c8-48c6-a7b2-8a69bd40fb44)

<div align="center">
    <h1>The Protein Ligand INteractions Dataset and Evaluation Resource</h1>
</div>

---

[![license](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](https://github.com/plinder-org/plinder/blob/master/LICENSE.txt)
[![publish](https://github.com/plinder-org/plinder/actions/workflows/main.yaml/badge.svg)](https://github.com/plinder-org/plinder/pkgs/container/plinder)
[![website](https://img.shields.io/badge/website-plinder-blue.svg)](https://www.plinder.sh/)
[![bioRxiv](https://img.shields.io/badge/bioRxiv-2024.07.17.603955-blue.svg)](https://www.biorxiv.org/content/10.1101/2024.07.17.603955)
[![docs](https://github.com/plinder-org/plinder/actions/workflows/docs.yaml/badge.svg)](https://plinder-org.github.io/plinder/)
[![coverage](https://github.com/plinder-org/plinder/raw/python-coverage-comment-action-data/badge.svg)](https://github.com/plinder-org/plinder/tree/python-coverage-comment-action-data)

![overview](https://github.com/user-attachments/assets/39d251b1-8114-4242-b9fc-e0cce900d22f)

# 📚 About

**PLINDER**, short for **p**rotein **l**igand **in**teractions **d**ataset and
**e**valuation **r**esource, is a comprehensive, annotated, high quality dataset and
resource for training and evaluation of protein-ligand docking algorithms:

- \> 400k PLI systems across > 11k SCOP domains and > 50k unique small molecules
- 750+ annotations for each system, including protein and ligand properties, quality,
  matched molecular series and more
- Automated curation pipeline to keep up with the PDB
- 14 PLI metrics and over 20 billion similarity scores
- Unbound \(_apo_\) and _predicted_ Alphafold2 structures linked to _holo_ systems
- _train-val-test_ splits and ability to tune splitting based on the learning task
- Robust evaluation harness to simplify and standard performance comparison between
  models.

The *PLINDER* project is a community effort, launched by the University of Basel,
SIB Swiss Institute of Bioinformatics, Proxima (formerly VantAI), NVIDIA, MIT CSAIL,
and will be regularly updated.

To accelerate community adoption, PLINDER will be used as the field’s new Protein-Ligand
interaction dataset standard as part of an exciting competition at the upcoming 2024
[Machine Learning in Structural Biology (MLSB)](https://mlsb.io#challenge) Workshop at NeurIPS, one of the field's premiere academic gatherings.
More details about the competition and other helpful practical tips can be found at our recent workshop repo:
[Moving Beyond Memorization](https://github.com/plinder-org/moving_beyond_memorisation).

### 👋 [Join the P(L)INDER user group Discord Server!](https://discord.gg/KgUdMn7TuS)


## 🔢 Plinder versions

We version the `plinder` dataset with two controls:

- `PLINDER_RELEASE`: the month stamp of the last RCSB sync
- `PLINDER_ITERATION`: value that enables iterative development within a release

We version the `plinder` application using an automated semantic
versioning scheme based on the `git` commit history.
The `plinder.data` package is responsible for generating a dataset
release and the `plinder.core` package makes it easy to interact
with the dataset.

#### 🐛🐛🐛 Known bugs:
- ~~Source dataset contains incorrect `entry_release_date` dates, please, use `query_index` to get correct dates patched.~~
- ~~Complexes containing nucleic acid receptors may [not be saved correctly](https://github.com/plinder-org/plinder/issues/61).~~
- ~~`ligand_binding_affinity` queries have been disabled due to a [bug found parsing BindingDB](https://github.com/plinder-org/plinder/issues/94)~~
All fixed in WIP — will take effect after dataset regeneration.

#### Changelog:

- WIP (Current — unreleased):
    - **Major backend refactor**: replaced OST, gemmi, plip, openbabel with biotite + peppr for data generation; removed 6 dependencies from ingest pipeline
    - **Nucleic acid support**: DNA/RNA chains now correctly included as receptor neighbors, mainchain/sidechain detection works for both protein and nucleic acids ([#61](https://github.com/plinder-org/plinder/issues/61))
    - **Custom CIF support**: parse Boltz/AF3/Chai outputs with bond order assignment from SMILES ([#117](https://github.com/plinder-org/plinder/issues/117))
    - **Stereochemistry**: CCD ideal 3D coordinates used as stereo ground truth; new `resolved_stereo_matches_template` flag validates resolved structure chirality against CCD template (handles partial resolution via MCS trimming)
    - **Interactions**: water bridge and metal bridge detection via peppr; halogen bond sidechain flag now computed (was hardcoded)
    - **Binding affinity**: fixed BindingDB matching — target sequence now validated against PDB SEQRES with 100% core identity, terminal tags/truncations tolerated ([#94](https://github.com/plinder-org/plinder/issues/94)); updated to BindingDB 2026-04
    - **Optional eval**: OpenStructure and posebusters moved to `pip install plinder[eval]`; base install is numpy 2 compatible; posebusters no longer runs during ingest
    - **PlinderSystem API**: new `receptor_structure` (biotite AtomArray) and `ligand_mols` (RDKit Mol) properties; OST properties (`receptor_entity`, `ligand_views`) kept for eval but require `plinder[eval]`
    - **Chain type support**: `Chain.from_cif_data` now assigns proper one-letter codes and chem_types for nucleotides (`RNA Linking`, `DNA Linking`); new `Residue.is_modified` property covers both protein PTMs and modified nucleotide bases
    - **Save utils**: receptor/ligand chain naming generalized (`PDB_RECEPTOR_CHAINS`); system saving works for protein, NA, and mixed complexes
    - **System definition**: unified `min_polymer_size=12` replaces separate `min_polymer_size`/`max_non_small_mol_ligand_length` — polymers ≥ 12 residues are receptor, shorter are ligands (threshold matches minimum MMseqs2/Foldseek search length); molecules with BIRD annotation are ligands irrespective of size; ligand chains no longer appear in both receptor and ligand parts of system IDs.
    - **System grouping**: pocket-based grouping (≥ 3 shared receptor residues on the same chain instance) for adjacent binding sites (e.g. orthosteric + allosteric, cofactor + substrate in same active site); artifacts attach only via 4 Å proximity
    - **Dead code removal**: removed unused OST-based functions, PDB string roundtrips, duplicate SMILES derivation paths, v1 template matching (consolidated to Rascal MCES `get_matched_template`)
    - **License**: changed from GPL-2.0 to Apache-2.0 (GPL was only required by PLIP, now removed)

- 2024-06/v2:
    - New systems added based on the 2024-06 RCSB sync
    - Updated system definition to be more stable and depend only on ligand distance rather than PLIP
    - Added annotations for crystal contacts
    - Improved ligand handling and saving to fix some bond order issues
    - Improved covalency detection and annotation to reference each bond explicitly
    - Added linked apo/pred structures to v2/links and v2/linked_structures
    - <del>Added binding affinity annotations from [BindingDB](https://bindingdb.org)</del> (see known bugs!)
    - Added statistics requirement and other changes in the split to enrich test set diversity

- 2024-04/v1: Version described in the preprint, with updated redundancy removal by protein pocket and ligand similarity.
- 2024-04/v0: Version used to re-train DiffDock in the paper, with redundancy removal based on \<pdbid\>\_\<ligand ccd codes\>

## 🏅 Gold standard benchmark sets

As part of *PLINDER* resource we provide train, validation and test splits that are
curated to minimize the information leakage based on protein-ligand interaction
similarity.
In addition, we have prioritized the systems that has a linked experimental `apo`
structure or matched molecular series to support realistic inference scenarios for hit
discovery and optimization.
Finally, a particular care is taken for test set that is further prioritized to contain
high quality structures to provide unambiguous ground-truths for performance
benchmarking.

![test_stratification](https://github.com/user-attachments/assets/5bb96534-f939-42b5-bf85-5ac3a71aa324)

Moreover, as we enticipate this resource to be used for benchmarking a wide range of methods, including those simultaneously predicting protein structure (aka. co-folding) or those generating novel ligand structures, we further stratified test (by novel ligand, pocket, protein or all) to cover a wide range of tasks.

# 👨💻 Getting Started

The *PLINDER* dataset is provided in two ways:

- You can either use the files from the dataset directly using your preferred tooling
  by downloading the data from the public
  [bucket](https://cloud.google.com/storage/docs/buckets),
- or you can utilize the dedicated `plinder` Python package for interfacing the data.


## Downloading the dataset

The dataset can be downloaded from the bucket with
[gsutil](https://cloud.google.com/storage/docs/gsutil_install).

```console
$ export PLINDER_RELEASE=2024-06 # Current release
$ export PLINDER_ITERATION=v2 # Current iteration
$ mkdir -p ~/.local/share/plinder/${PLINDER_RELEASE}/${PLINDER_ITERATION}/
$ gsutil -m cp -r "gs://plinder/${PLINDER_RELEASE}/${PLINDER_ITERATION}/*" ~/.local/share/plinder/${PLINDER_RELEASE}/${PLINDER_ITERATION}/
```
For details on the sub-directories, see [Documentation](https://plinder-org.github.io/plinder/tutorial/dataset.html).

## Installing the Python package

`plinder` is available on *PyPI*.

```
pip install plinder
```

For evaluation scoring (lDDT, RMSD via OpenStructure):

```
pip install plinder[eval]
```

## License
Data curated by PLINDER are made available under the Apache License 2.0.
All data curated by BindingDB staff are provided under the Creative Commons Attribution 4.0 License. Data imported from ChEMBL are provided under their Creative Commons Attribution-Share Alike 4.0 Unported License.


# 📝 Documentation

A more detailed description is available on the
[documentation website](https://plinder-org.github.io/plinder/).

# 📃 Citation

Durairaj, Janani, Yusuf Adeshina, Zhonglin Cao, Xuejin Zhang, Vladas Oleinikovas, Thomas Duignan, Zachary McClure, Xavier Robin, Gabriel Studer, Daniel Kovtun, Emanuele Rossi, Guoqing Zhou, Srimukh Prasad Veccham, Clemens Isert, Yuxing Peng, Prabindh Sundareson, Mehmet Akdel, Gabriele Corso, Hannes Stärk, Gerardo Tauriello, Zachary Wayne Carpenter, Michael M. Bronstein, Emine Kucukbenli, Torsten Schwede, Luca Naef. 2024. “PLINDER: The Protein-Ligand Interactions Dataset and Evaluation Resource.”
[bioRxiv](https://doi.org/10.1101/2024.07.17.603955)
[ICML'24 ML4LMS](https://openreview.net/forum?id=7UvbaTrNbP)

Please see the [citation file](CITATION.cff) for details.

![plinder_banner](https://github.com/user-attachments/assets/43d129f2-3bb6-4903-81fa-182c351c64b6)
