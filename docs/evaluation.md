---
sd_hide_title: true
---

# Evaluation

Compare predicted ligand poses and protein interfaces with PLINDER references
using OpenStructure. PoseBusters checks the chemical and physical plausibility
of predicted ligands.

## Installation

Install OpenStructure 2.12.0 or newer from Bioconda and the Python evaluation extra:

```bash
conda install -c conda-forge -c bioconda 'openstructure>=2.12.0'
pip install 'plinder[eval]'
```

The `ost` executable must be on your `PATH`. PLINDER calls its
`compare-ligand-structures` and `compare-structures` actions directly.
The evaluation tests use OpenStructure 2.12.0 and PoseBusters 0.6.5.

## Arrange your predictions

Put each complete predicted mmCIF in a subdirectory named after its reference:

```text
predictions/
├── 1avd__1__1.A__1.C/
│   ├── model_1.cif
│   └── model_2.cif
└── 1avd/
    └── model_1.cif
```

Use a ligand-system ID or protein-interface ID for a specific reference.
A four-character PDB ID compares each prediction with all matching systems
or interfaces from that entry. Each file is a separate prediction.
Both `.cif` and `.mmcif` files are accepted, including gzip-compressed files;
the first coordinate model in each file is used.

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
comparisons in sequence; PoseBusters does not start another process pool.

## Ligand selection and bonds

Ligand evaluation uses **proper reference ligands** by default: those marked
`ligand_is_proper` in the annotation table. Add `--include-all-ligands` to
also evaluate reference ions and artifacts. OST assigns predicted ligands to
the selected references. Predictions are not discarded for being distant from
the receptor.

PLINDER prepares ligand SDFs from the complete prediction. Existing bonds and
CCD chemistry are used where available. For an unfamiliar component name such
as `LIG`, supply its CCD code or SMILES:

```bash
plinder_eval predictions --output-dir evaluation --mode ligands \
    --ligand-ccd-codes '{"LIG": "ATP"}'
```

Alternatively use `--ligand-smiles '{"LIG": "CCO"}'`.
For SMILES, the heavy-atom order must match the CIF atom order. CCD templates
use atom names or graph matching. These mappings apply to every prediction
in the run. Supply only one chemistry source per component.

Protein and DNA/RNA polymers are treated as receptors. For a predicted peptide
ligand, use `--ligand-chain B`, where `B` is its label asym ID; repeat the
option for additional polymer ligand chains. DNA/RNA inputs need their mmCIF
polymer sequence metadata. Covalently connected ligand chains stay together;
metal coordination does not combine separate ligands.

PoseBusters is enabled for ligand evaluation. Use `--no-posebusters` to skip it.
The receptor retains its residue information in memory, without a PDB conversion.

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
original label asym IDs; a covalent multi-chain group uses its lexically first ID.

Interface metrics include `lddt`, `ilddt`, `qs_global`, `qs_best`,
`dockq`, `dockq_ave_full`, and `dockq_wave_full`. The full DockQ aggregates
include zero penalties for unmapped reference interfaces.

Unmatched reference ligands remain in the table, with a reason for each
unassigned metric. Their scores are unavailable, not invented zeros.
`status="unassigned"` means neither ligand metric found an assignment.
Include these rows as failures when calculating a pose-success fraction.
`status="error"` indicates an execution problem; inspect `failures.tsv`
before interpreting results. These errors are distinct from poor predictions.

The command returns a nonzero exit code if errors were recorded; unassigned
ligands alone do not cause a nonzero exit code. Each run recalculates results
and replaces the summary tables. Old detail files are never used as results.

For additional OST settings, pass a quoted argument string with
`--ligand-options='--flag value'` or `--interface-options='--flag value'`.
The Python function accepts the corresponding lists of arguments.
