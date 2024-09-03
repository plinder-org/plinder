---
sd_hide_title: true
html_theme.sidebar_secondary.remove: true
---

# PLINDER documentation


![plinder](/static/assets/general/plinder_logo.png){w=40em align=center}

**PLINDER**, short for **p**rotein **l**igand **in**teractions **d**ataset and
**e**valuation **r**esource, is a comprehensive, annotated, high quality dataset and
resource for training and evaluation of protein-ligand docking algorithms:

- \> 400k PLI systems across > 11k SCOP domains and > 50k unique small molecules
- 500+ annotations for each system, including protein and ligand properties, quality,
  matched molecular series and more
- Automated curation pipeline to keep up with the PDB
- 14 PLI metrics and over 20 billion similarity scores
- Unbound \(_apo_\) and _predicted_ Alphafold2 structures linked to _holo_ systems
- _train-val-test_ splits and ability to tune splitting based on the learning task
- Robust evaluation harness to simplify and standard performance comparison between models.


::::::{grid} 1 1 2 2

:::::{grid-item-card}
:link: tutorial/dataset
:link-type: doc

::::{grid} 2

:::{grid-item}
:columns: 3
:class: main-button
<svg class="main-button-icon" xmlns="http://www.w3.org/2000/svg" viewBox="0 0 384 512"><!--!Font Awesome Free 6.6.0 by @fontawesome - https://fontawesome.com License - https://fontawesome.com/license/free Copyright 2024 Fonticons, Inc.--><path d="M0 64C0 28.7 28.7 0 64 0L224 0l0 128c0 17.7 14.3 32 32 32l128 0 0 288c0 35.3-28.7 64-64 64L64 512c-35.3 0-64-28.7-64-64L0 64zm384 64l-128 0L256 0 384 128z"/></svg>
:::

:::{grid-item}
:columns: 9
**Dataset access**

Access the PLI systems and their annotations directly via the files
:::
::::
:::::


:::::{grid-item-card}
:link: tutorial/api
:link-type: doc

::::{grid} 2

:::{grid-item}
:columns: 3
:class: main-button
<svg class="main-button-icon" xmlns="http://www.w3.org/2000/svg" viewBox="0 0 640 512"><!--!Font Awesome Free 6.6.0 by @fontawesome - https://fontawesome.com License - https://fontawesome.com/license/free Copyright 2024 Fonticons, Inc.--><path d="M392.8 1.2c-17-4.9-34.7 5-39.6 22l-128 448c-4.9 17 5 34.7 22 39.6s34.7-5 39.6-22l128-448c4.9-17-5-34.7-22-39.6zm80.6 120.1c-12.5 12.5-12.5 32.8 0 45.3L562.7 256l-89.4 89.4c-12.5 12.5-12.5 32.8 0 45.3s32.8 12.5 45.3 0l112-112c12.5-12.5 12.5-32.8 0-45.3l-112-112c-12.5-12.5-32.8-12.5-45.3 0zm-306.7 0c-12.5-12.5-32.8-12.5-45.3 0l-112 112c-12.5 12.5-12.5 32.8 0 45.3l112 112c12.5 12.5 32.8 12.5 45.3 0s12.5-32.8 0-45.3L77.3 256l89.4-89.4c12.5-12.5 12.5-32.8 0-45.3z"/></svg>
:::

:::{grid-item}
:columns: 9
**Python API**

Use the dedicated Python package to explore the data
:::
::::
:::::

::::::

% TODO: re-add `contribution/index`

:::{toctree}
:maxdepth: 1
:hidden:

tutorial/index
dataset
api/index
evaluation
examples/index
contribution/index
citation
:::
