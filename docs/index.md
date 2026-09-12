---
sd_hide_title: true
html_theme.sidebar_secondary.remove: true
---

# PLINDER documentation

![plinder](/static/assets/general/plinder_logo.png){w=40em align=center}

**PLINDER** is the **Protein & Ligand INteraction Dataset and Evaluation
Resource**: a comprehensive, annotated, high-quality resource for training and
evaluating protein-ligand and protein-protein structure models.

- \> 400k PLI systems across > 11k SCOP domains and > 50k unique small molecules
- Ligand and protein-interface annotations plus compact entry, chain, and representative tables
- Automated curation pipeline to keep up with the PDB
- Reusable ligand, pocket, and protein-interface similarities and representative assignments
- Ranked deposited apo protein chains linked to _holo_ ligands
- Python APIs for release queries, coordinate reconstruction, and custom scoring

::::::{grid} 1 2 3 3

:::::{grid-item-card} Getting started
:link: tutorial/dataset
:link-type: doc
:class-card: home-card

<div class="home-card-icon" aria-hidden="true"><i class="fa-solid fa-compass"></i></div>

Install PLINDER, select a release, and continue to the runnable guides.
:::::

:::::{grid-item-card} Data access
:link: examples/2_query_filter_index
:link-type: doc
:class-card: home-card

<div class="home-card-icon" aria-hidden="true"><i class="fa-solid fa-database"></i></div>

Download release files and query columns across the ligand, interface,
entry, chain, and representative tables.
:::::

:::::{grid-item-card} Structures
:link: examples/3_access_system_files
:link-type: doc
:class-card: home-card

<div class="home-card-icon" aria-hidden="true"><i class="fa-solid fa-cubes"></i></div>

Work with ligand systems, linked apo chains, protein interfaces, atom masks,
and self-contained mmCIF files.
:::::

:::::{grid-item-card} Similarity and scoring
:link: examples/index
:link-type: doc
:class-card: home-card

<div class="home-card-icon" aria-hidden="true"><i class="fa-solid fa-code-compare"></i></div>

Use pairwise scores and representative covers, or compare custom
mmCIF and FASTA inputs with PLINDER.
:::::

:::::{grid-item-card} Data reference
:link: dataset
:link-type: doc
:class-card: home-card

<div class="home-card-icon" aria-hidden="true"><i class="fa-solid fa-book-open"></i></div>

Inspect release files, table relationships, column definitions, and
coordinate provenance.
:::::

::::::

:::{toctree}
:maxdepth: 1
:hidden:

tutorial/index
examples/index
dataset
evaluation
api/index
contribution/index
citation
:::
