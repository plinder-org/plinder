---
sd_hide_title: true
---

# Evaluation

Compare predicted ligand poses and protein interfaces with PLINDER references
using OpenStructure. PoseBusters checks the chemical and physical plausibility
of predicted ligands.

For a runnable walkthrough, see {doc}`examples/evaluation`.

## Installation

Install OpenStructure 2.12.0 or newer from Bioconda and the Python evaluation extra:

```bash
conda install -c conda-forge -c bioconda 'openstructure>=2.12.0'
pip install 'plinder[eval]'
```

The `ost` executable must be on your `PATH`. PLINDER calls its
`compare-ligand-structures` and `compare-structures` actions directly.
The evaluation tests use OpenStructure 2.12.0 and PoseBusters 0.6.5.

(arrange-your-predictions)=
## Arrange your predictions

Search and evaluation accept the same CSV, TSV, Parquet or pandas input table:

| input_id | structure_path | ligand_path | reference_id |
| --- | --- | --- | --- |
| model_1 | models/receptor.pdb | models/pose.sdf | 8c3u__1__1.A__1.C |
| model_2 | models/complex.cif | | 8c3u__1__1.A__1.C |
| model_3 | models/interface.pdb | | 2e31__1__1.A--1.B |

`input_id` identifies a prediction. For multiple ligand SDFs, repeat its row
with a different `ligand_path`, or point `ligand_path` to a folder of SDFs.
An empty ligand field selects a complete structure; an empty ligand folder
represents zero predicted poses. Paths are relative to the table's directory;
DataFrames use the working directory. Evaluation requires `reference_id`,
which can be a PLINDER system/interface ID or a PDB ID. Search uses the remaining
columns to find matches across PLINDER.

Pass the table to `evaluate_predictions("inputs.tsv", output_dir="evaluation")`
or `plinder_eval inputs.tsv --output-dir evaluation`. The search equivalent is
`plinder.core.scores.search("inputs.tsv", output_dir="search_results")`;
see {doc}`examples/custom_scoring` for a walkthrough.

### Folder input

Put each prediction in a subdirectory named after its reference:

```text
predictions/
├── 8c3u__1__1.A__1.C/
│   ├── model_1.pdb
│   ├── model_1.sdf
│   ├── model_2.cif
│   └── model_2.ligands/
│       ├── ligand_A.sdf
│       └── ligand_B.sdf
└── 8c3u/
    └── model_1.cif
```

Use a ligand-system ID or protein-interface ID for a specific reference.
A four-character PDB ID compares each prediction with all matching systems
or interfaces from that entry. Each file is a separate prediction.
Both `.cif` and `.mmcif` files are accepted, including gzip-compressed files;
the first coordinate model in each file is used.
For protein interfaces, `.pdb` and `.pdb.gz` predictions can also be used directly.

For ligand evaluation, supply a predicted receptor PDB/mmCIF with a matching
`model_1.sdf`, or a `model_1.ligands/` folder containing one SDF per ligand.
The receptor and ligand coordinates must describe the same predicted complex.
Each SDF supplies its ligand's chemistry and coordinates. An empty
`model_1.ligands/` folder represents a prediction with zero ligand poses.
Choose either the single SDF or the ligand folder for each receptor.

A complete predicted mmCIF can also be used on its own; PLINDER extracts its
receptor and ligands for evaluation.

## Run evaluation

```bash
plinder_eval predictions --output-dir evaluation --mode both --num-workers 8
```

Use `--mode ligands` or `--mode interfaces` to evaluate only one type.
References come from your configured PLINDER release. To use a local release,
add `--data-dir /path/to/release`. Reference coordinates are reconstructed
from the release's recorded source structures as needed.

The equivalent Python function is
`plinder.eval.evaluate_predictions(predictions, output_dir=..., mode="both", num_workers=8)`.
It returns a dictionary of DataFrames as well as saving them to the output folder.
In a standalone Python script, put the call inside
`if __name__ == "__main__":`.

`num_workers` limits concurrent prediction processes. Each worker runs its
comparisons and PoseBusters checks in sequence.

## Ligand selection and bonds

Ligand evaluation uses **proper reference ligands** by default: those marked
`ligand_is_proper` in the annotation table. Add `--include-all-ligands` to
also evaluate reference ions and artifacts. OST assigns predicted ligands to
the selected references, including predicted ligands far from the receptor.

For receptor + SDF inputs, chemistry comes from each prediction's own SDF files.
Use this layout when different predictions use `LIG` for different molecules.
Results identify these ligands by their original SDF filenames.

For complete-mmCIF inputs, PLINDER prepares ligand SDFs from the prediction.
Existing bonds and CCD chemistry are used where available. For an unfamiliar
component name such as `LIG`, supply its CCD code or SMILES:

```bash
plinder_eval predictions --output-dir evaluation --mode ligands \
    --ligand-ccd-codes '{"LIG": "ATP"}'
```

Alternatively use `--ligand-smiles '{"LIG": "CCO"}'`.
For SMILES, the heavy-atom order must match the CIF atom order. CCD templates
use atom names or graph matching. These mappings apply to every prediction
in the run that uses complete-mmCIF extraction. Use a shared mapping only when
the component name refers to the same molecule across those predictions.
Supply only one chemistry source per component.

Protein and DNA/RNA polymers are treated as receptors. For a predicted peptide
ligand, use `--ligand-chain B`, where `B` is its label asym ID; repeat the
option for additional polymer ligand chains. DNA/RNA inputs need their mmCIF
polymer sequence metadata. Covalently connected ligand chains stay together;
ligands connected only by metal coordination are treated as separate molecules.

PoseBusters is enabled for ligand evaluation. Use `--no-posebusters` to skip it.
It checks each ligand against its paired predicted receptor, or the receptor
extracted from the same complete mmCIF.
The receptor's residue information is preserved in an RDKit molecule.

Covalently bound ligands need care when interpreting these checks: separating
the ligand and receptor can make their covalent attachment appear as a
nonbonded clash. Such checks remain `False` in the results; inspect the attachment
before interpreting it as a bad pose.

## Results

| File | Contents |
| --- | --- |
| `ligands.parquet` | One row per selected reference ligand and prediction |
| `interfaces.parquet` | One row per reference interface and prediction |
| `posebusters.parquet` | Plausibility checks per candidate predicted ligand |
| `failures.tsv` | Reference-loading, preparation, or tool errors |
| `details/` | Native OST results, PoseBusters CSVs, logs, and prepared structures |

The Python return keys are `ligands`, `interfaces`, `posebusters`, and `failures`.
`prediction` is the input path relative to the predictions folder;
`system_id` and `ligand_id` identify the reference.

Ligand metrics include `rmsd` (binding-site-superposed, symmetry-corrected RMSD),
`bb_rmsd`, `lddt_lp`, `lddt_pli`, and the corresponding ligand coverage values.
OST can choose different assignments for RMSD and lDDT-PLI. Their predicted
ligand IDs are therefore stored separately as `rmsd_model_ligand` and
`lddt_pli_model_ligand`.

PoseBusters includes all candidate predicted ligands, including unmatched ones.
To attach its checks to a reference-ligand result, join on `prediction` and
the model-ligand ID for the metric you are analysing. The predicted IDs are
SDF filenames for supplied poses, or label asym IDs for complete-mmCIF inputs;
a covalent multi-chain group uses its lexically first asym ID.

Interface metrics include `lddt`, `ilddt`, `qs_global`, `qs_best`,
`dockq`, `dockq_ave_full`, and `dockq_wave_full`. The full DockQ aggregates
include zero penalties for unmapped reference interfaces.

Unmatched reference ligands have empty scores and a reason for each
unassigned metric. `status="unassigned"` means both ligand metrics are unassigned.
Include these rows as failures when calculating a pose-success fraction.
`status="error"` indicates an execution problem; inspect `failures.tsv`
and resolve it before interpreting results.
`status="success"` means the comparison completed.

OpenStructure 2.12 can reject modified receptor residues when its residue code
disagrees with the sequence in the mmCIF, for example `CSK` in PDB entry `4agi`.
This appears as a sequence-mismatch error in `failures.tsv`.
Keep track of these cases when reporting how many predictions were evaluated.

The command returns exit code 0 when comparisons complete, including those with
unassigned ligands, and a nonzero code for execution errors. Each run recalculates
results and replaces the summary tables.

For additional OST settings, pass a quoted argument string with
`--ligand-options='--flag value'` or `--interface-options='--flag value'`.
The Python function accepts the corresponding lists of arguments.
