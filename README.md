![plinder](https://github.com/user-attachments/assets/05088c51-36c8-48c6-a7b2-8a69bd40fb44)

<div align="center">
    <h1>Protein &amp; Ligand INteraction Dataset and Evaluation Resource</h1>
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

**PLINDER** is the **Protein & Ligand INteraction Dataset and Evaluation
Resource**: a comprehensive, annotated, high-quality resource for training and
evaluating protein-ligand and protein-protein structure models.

- \> 400k PLI systems across > 11k SCOP domains and > 50k unique small molecules
- Ligand-level annotations plus compact entry, chain, interface, and representative tables
- Automated curation pipeline to keep up with the PDB
- Reusable ligand, pocket, and protein-interface similarities and cluster assignments
- Deposited apo protein chains linked to compatible _holo_ systems
- Directed-cover assignments for choosing diverse training representatives
- Robust evaluation harness to simplify and standard performance comparison between
  models.

The *PLINDER* project is a community effort, launched by the University of Basel,
SIB Swiss Institute of Bioinformatics, Proxima (formerly VantAI), NVIDIA, MIT CSAIL,
and will be regularly updated.

PLINDER set a new standard for the Protein-Ligand interaction datasets. It was first introduced as part of the 2024 Machine Learning in Structural Biology (MLSB) [Workshop challenge](https://www.mlsb.io/index_2024.html#challenge) at NeurIPS, one of the field's premiere academic gatherings.
More details about the competition and other helpful practical tips can be found at our recent workshop repo:
[Moving Beyond Memorization](https://github.com/plinder-org/moving_beyond_memorisation).

### 👋 [Join the P(L)INDER user group Discord Server!](https://discord.gg/KgUdMn7TuS)


## 🔢 Plinder versions

We version the `plinder` dataset with two controls:

- `PLINDER_RELEASE`: the month stamp of the last RCSB sync
- `PLINDER_RELEASE_NUMBER`: numbered release within that ingest month

We version the `plinder` application using an automated semantic
versioning scheme based on the `git` commit history.
The `plinder.data` package is responsible for generating a dataset
release and the `plinder.core` package makes it easy to interact
with the dataset.

#### Changelog:

- WIP (Current — unreleased):
    - **Major backend refactor**: replaced OST, gemmi, plip, openbabel with biotite + peppr for data generation; removed 6 dependencies from ingest pipeline
    - **Nucleic acid support**: DNA/RNA chains now correctly included as receptor neighbors, mainchain/sidechain detection works for both protein and nucleic acids ([#61](https://github.com/plinder-org/plinder/issues/61))
    - **Custom CIF support**: new `Entry.from_custom_cif_file` for structure-prediction outputs (Boltz, AlphaFold3, Chai-1) that ship CIFs without `_chem_comp_bond` ([#117](https://github.com/plinder-org/plinder/issues/117)). Bond orders come from `ligand_smiles_dict` via positional atom-order match (the convention these tools follow); element/count mismatches raise with the offending position, `force_substructure_match=True` opts into substructure matching when atom order isn't preserved. User SMILES win over CCD for both `smiles` and `resolved_stereo_matches_template` — closes a silent gap where biotite's `LIG` placeholder would pass any 3D conformer. Input CIFs are never mutated; optional `save_fixed_cif` persists the enriched copy.
    - **Stricter CIF ingest**: H/D filtered consistently (biotite's `filter_heavy`); multi-model CIFs warn and use model 1; multi-instance custom comp_ids must share heavy-atom naming (since `_chem_comp_bond` is comp_id-keyed); silent `connect_via_residue_names` and half-sanitized substructure fallbacks replaced with `ValueError` so corrupt inputs fail loudly.
    - **Stereochemistry**: CCD ideal 3D coordinates used as stereo ground truth; new `resolved_stereo_matches_template` flag validates resolved structure chirality against CCD template (handles partial resolution via MCS trimming)
    - **Interactions**: water bridge and metal bridge detection via peppr; halogen bond sidechain flag now computed (was hardcoded)
    - **Binding affinity**: fixed BindingDB matching — target sequence now validated against PDB SEQRES with 100% core identity, terminal tags/truncations tolerated ([#94](https://github.com/plinder-org/plinder/issues/94)); updated code to get the latest BindingDB release
    - **Optional eval**: `pip install plinder[eval]` adds PoseBusters ligand validation; OpenStructure-backed metrics use the Conda-only `openstructure` package; PoseBusters no longer runs during ingest
    - **PlinderSystem API**: new `receptor_structure` and `ligand_structures` (Biotite AtomArray) plus `ligand_mols` (RDKit Mol) properties; OpenStructure is confined to the evaluation implementation
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

## 🏅 Preprint benchmark sets

The historical `2024-04/v1` preprint release provides train, validation, and test splits that are
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

After installing the package, download the index tables and cluster assignments
for a release with:

```console
$ plinder_download --release 2026-07 --release-number 1
```

The command offers the larger ligand, alignment, score, export, and search
database groups separately. Missing optional artifacts are fetched when an API
call needs them. Files can also be copied directly from the public bucket with
[`gsutil`](https://cloud.google.com/storage/docs/gsutil_install).
For details on the release paths, see [Documentation](https://plinder-org.github.io/plinder/tutorial/dataset.html).

## Installing the Python package

`plinder` is available on *PyPI*.

```
pip install plinder
```

For PoseBusters ligand validation:

```
pip install plinder[eval]
```

OpenStructure is not published on PyPI. Evaluation metrics backed by
OpenStructure (including lDDT and RMSD) use its command-line actions and require
OpenStructure 2.12.0 or newer; the repository's `environment.yml` installs it
from Bioconda. See the [evaluation guide](docs/evaluation.md) for evaluating
folders of ligand or protein-interface predictions.

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
