![plinder](https://github.com/user-attachments/assets/05088c51-36c8-48c6-a7b2-8a69bd40fb44)

# PLINDER

Protein & Ligand INteraction Dataset and Evaluation Resource

[![license](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](LICENSE.txt)
[![publish](https://github.com/plinder-org/plinder/actions/workflows/main.yaml/badge.svg)](https://pypi.org/project/plinder/)
[![docs](https://github.com/plinder-org/plinder/actions/workflows/docs.yaml/badge.svg)](https://plinder-org.github.io/plinder/)
[![bioRxiv](https://img.shields.io/badge/bioRxiv-2024.07.17.603955-blue.svg)](https://doi.org/10.1101/2024.07.17.603955)

PLINDER brings together protein and nucleic-acid monomers, protein–ligand systems,
protein–protein interfaces, protein–nucleic-acid complexes, and their
similarities for training and evaluating structure prediction and cofolding models.

The September 2026 release contains:

<!-- Regenerate: python scripts/update_readme_counts.py /path/to/release --readme README.md -->
<!-- release-counts:start -->
| Release contents | Count |
| --- | ---: |
| Ligand systems | 611,930 |
| Proper ligands | 803,190 |
| Protein–protein interfaces | 3,811,623 |
| Protein-monomer assemblies | 161,654 |
| Nucleic-acid-monomer assemblies | 2,208 |
| Protein–nucleic-acid assemblies | 16,766 |

Nucleic-acid monomers: 600 DNA, 1,587 RNA, and 21 DNA/RNA hybrids.
<!-- release-counts:end -->

## Start with the Python API

```bash
pip install plinder
```

Files are downloaded as you use them and cached locally. Start with a small
query; downloading the full release is optional.

```python
from plinder.core import query_table

ligands = query_table(
    "annotation",
    columns=["system_id", "ligand_id", "ligand_unique_ccd_code", "entry_resolution"],
    filters=[("entry_pdb_id", "==", "4agi"), ("ligand_is_proper", "==", True)],
)
print(ligands)
```

For protein–protein complexes, query `"interface_annotations"`. Use
`PlinderSystem` or `PlinderInterface` to access sequences, atom arrays,
ligand molecules, and reconstructed mmCIF files.

## What can I do with PLINDER?

| Task | Runnable guide |
| --- | --- |
| Select ligands, interfaces, chains, and experimental annotations | [Query the whole annotated PDB](https://plinder-org.github.io/plinder/examples/2_query_filter_index.html) |
| Access protein–ligand systems and linked apo structures as CIF, SDF, FASTA, SMILES, and arrays | [Protein–ligand systems](https://plinder-org.github.io/plinder/examples/3_access_system_files.html) |
| Work with coordinate masks and ligand atom correspondences | [Atom arrays and mappings](https://plinder-org.github.io/plinder/examples/4_align_mask_crop.html) |
| Access protein–protein interfaces and their apo chains | [Protein–protein interfaces](https://plinder-org.github.io/plinder/examples/protein_interfaces.html) |
| Query ligand, pocket, chain, whole-interface, and half-interface scores | [Similarity scores](https://plinder-org.github.io/plinder/examples/similarity_and_representatives.html) |
| Select representatives and examine clusters | [Representative covers and components](https://plinder-org.github.io/plinder/examples/representatives_and_leakage.html) |
| Search with your own protein or complex sequences and structures | [Custom searches](https://plinder-org.github.io/plinder/examples/custom_scoring.html) |
| Evaluate ligand poses and protein-interface predictions | [Evaluation](https://plinder-org.github.io/plinder/examples/evaluation.html) |
| Download the full release locally | [Downloads](https://plinder-org.github.io/plinder/examples/1_download.html) |

The [dataset reference](https://plinder-org.github.io/plinder/dataset.html)
describes release files, table relationships, and columns. The
[Python API reference](https://plinder-org.github.io/plinder/api/index.html)
lists individual functions and classes.

## Search and evaluate your models

`plinder.core.scores.search()` accepts protein FASTA files, structure files or
folders, and tables pairing receptors with ligand SDFs. It returns chain hits
and the applicable pocket, ligand, or interface scores. CIF/PDB searches use
MMseqs2 and Foldseek; FASTA queries use MMseqs2.

`plinder.eval.evaluate_predictions()` compares predictions with release
references and returns ligand/interface metrics and PoseBusters checks.
Search and evaluation use the same structure-table convention.

Custom searches require the corresponding MMseqs2/Foldseek executables.
For evaluation, install `plinder[eval]` and OpenStructure 2.12.0 or newer;
see the [installation instructions](https://plinder-org.github.io/plinder/evaluation.html#installation).

## Releases and local files

The current dataset release is `2026-09`. Package versions are separate from
dataset releases. See [CHANGELOG.md](CHANGELOG.md) for the release history.

Public files are available [here](https://cameo3d.org/plinder/PLINDER-2026-09/).
Set `PLINDER_MOUNT` before importing PLINDER to choose a cache location.
Use `plinder_download` to prepare a local copy and `PLINDER_OFFLINE=true`
to use cached files without network access.

## Community and examples

PLINDER is a community effort launched by the University of Basel,
SIB Swiss Institute of Bioinformatics, Proxima (formerly VantAI), NVIDIA,
and MIT CSAIL.

Join the [P(L)INDER Discord](https://discord.gg/KgUdMn7TuS), or explore the
[Moving Beyond Memorization workshop](https://github.com/plinder-org/moving_beyond_memorisation)
for the original 2024 tutorials on dataset selection, similarity, and evaluation.
The guides linked above use the current API and release.

## License

PLINDER code and curated data are available under the Apache License 2.0.
BindingDB-curated data use CC BY 4.0; data imported from ChEMBL use
CC BY-SA 4.0.

## Citation

Durairaj, Janani, Yusuf Adeshina, Zhonglin Cao, Xuejin Zhang, Vladas Oleinikovas, Thomas Duignan, Zachary McClure, Xavier Robin, Gabriel Studer, Daniel Kovtun, Emanuele Rossi, Guoqing Zhou, Srimukh Prasad Veccham, Clemens Isert, Yuxing Peng, Prabindh Sundareson, Mehmet Akdel, Gabriele Corso, Hannes Stärk, Gerardo Tauriello, Zachary Wayne Carpenter, Michael M. Bronstein, Emine Kucukbenli, Torsten Schwede, Luca Naef. 2024. “PLINDER: The Protein-Ligand Interactions Dataset and Evaluation Resource.”
[bioRxiv](https://doi.org/10.1101/2024.07.17.603955)
[ICML'24 ML4LMS](https://openreview.net/forum?id=7UvbaTrNbP)

Please see the [citation file](CITATION.cff) for details.
