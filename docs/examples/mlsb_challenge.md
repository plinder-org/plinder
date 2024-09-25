# MLSB/P(L)INDER Challenge: Overview and details

## Objectives

**Benchmarking:** Establish a standardized evaluation framework for protein-protein and protein-ligand structure prediction tasks. <br>
**Progress Tracking:** Track advancements in the field over time through a leaderboard that provides an objective basis for comparing different methods. <br>

## Key Dates

**First Announcement:** September 4th, 2024. <br>
**Detailed Instructions Release:** September 24th, 2024, at the [training workshop](https://unibas.zoom.us/meeting/register/u5EkcOusrzsiHtLDImocB5PM0RLK9vC-g-kW#/registration).<br>
**Leaderboard Opens:** October 9th, 2024 (following acceptance notifications for MLSB).<br>
**Leaderboard Closes:** November 9th, 2024.<br>
**Winner Notification:** Wednesday, November 27th, 2024.<br>
**MLSB Workshop:** December 14-15th, 2024.<br>

(mlsb-rules-target)=

## Rules for valid model training

- Participants **MUST** use the sequences and SMILES in the provided train and validation sets from PINDER or PLINDER. In order to ensure no leakage, external data augmentation is not allowed.
- If starting structures/conformations need to be generated for the model, then this can only be done from the training and validation sequences and SMILES. Note that this is only the case for train & validation - no external folding methods or starting structures are allowed for the test set under any circumstance!. Only the predicted structures/conformers themselves may be used in this way, the embeddings or models used to generate such predictions may not. E.g. it is not valid to “distill” a method that was not trained on PLINDER/PINDER
- The PINDER and PLINDER datasets should be used independently; combining the sets is considered augmentation and is not allowed.
- For inference, only the inputs provided in the evaluation sets may be used: canonical sequences, structures and MSAs; no alternate templates or sequences are permitted. The inputs that will be used by assessors for each challenge track is as follows:
  - PLINDER: (SMILES, monomer protein structure, monomer FASTA, monomer MSA)
  - PINDER: (monomer protein structure 1, monomer protein structure 2, FASTA 1, FASTA 2, MSA 1, MSA 2)
- Model selection must be performed exclusively on the validation set designed for this purpose within the PINDER and PLINDER datasets.
- Methods relying on any model derivatives or embeddings trained on structures outside the PINDER/PLINDER training set are not permitted (e.g., ESM2, MSA: ✅; ESM3/ESMFold/SAProt/UniMol: ❌).
- For instruction on how to load training and validation data, check the links below:
  - [PLINDER](https://plinder-org.github.io/plinder/examples/index.html)
  - [PINDER](https://pinder-org.github.io/pinder/pinder-mlsb.html#accessing-and-loading-data-for-training)

## Rules for valid inference pipeline

- Using additional template structures or sequences is NOT allowed
- No binding site information can be used
- For each system, you are allowed to submit at most 1 prediction! If your method produces multiple samples, you must rank/score the predictions and only supply a single prediction as the top-ranking prediction that will be used in the leaderboard.
- Model inference should run in under 10 minutes per system on a GPU like T4, A10G
- The final predictions must be in:
  - PINDER track: CIF/PDB file format and contain two chains: Receptor (chain R) and Ligand (chain L)
  - PLINDER track: CIF/PDB file format for Receptor (protein) and SDF file format for Ligand (small molecule)
- Systems without a valid prediction will be penalized.

## Rules for valid submission

Submission system will use Hugging Face Spaces. To qualify for submission, each team must:

- Provide an MLSB submission ID or a link to a preprint/paper describing their methodology. This publication does not have to specifically report training or evaluation on the P(L)INDER dataset. Previously published methods, such as DiffDock, only need to link their existing paper. Note that entry into this competition does not equate to an MLSB workshop paper submission.
- Create a copy of the provided [inference template](https://huggingface.co/spaces/MLSB/plinder_inference_template/blob/main/inference_app.py).
  - Go to the top right corner of the page and click on the drop-down menu (vertical ellipsis) right next to the “Community”, then select “Duplicate this space”.
- Change files in the newly create space to reflect the peculiarities of your model
  - Edit `requirements.txt` to capture all dependencies.
  - Include a `inference_app.py` file. This contains a `predict` function that should be modified to reflect the specifics of inference using their model.
  - Include a `train.py` file to ensure that training and model selection use only the PINDER/PLINDER datasets and to clearly show any additional hyperparameters used.
  - Provide a LICENSE file that allows for reuse, derivative works, and distribution of the provided software and weights (e.g., MIT or Apache2 license).
  - Modify the Dockerfile as appropriate (including selecting the right base image)
- Submit to the leaderboard via the [designated form](https://huggingface.co/spaces/MLSB/leaderboard2024).
  - On submission page, add reference to the newly created space in the format username/space (e.g mlsb/alphafold3)

## Metrics

- Primary Ranking Metric:
  - PLINDER: lDDT-PLI
  - PINDER: DockQ

Other metrics computed by PINDER/PLINDER will be displayed on the leaderboard but will not influence the ranking.

## Evaluation Datasets

Although the exact composition of the eval set will be shared at a future date, below we provide an overview of the dataset and what to expect

- Two leaderboards, one for each of PINDER and PLINDER, will be created using a single evaluation set for each.
- Evaluation sets will be subsets of 150-200 structures from the current PINDER and PLINDER test splits (subsets to enable reasonable eval runtime).
- Each evaluation sample will contain a predefined input/output to ensure performance assessment is model-dependent, not input-dependent.
- The focus will be exclusively on flexible docking/co-folding, with a single canonical structure per protein, sampled from apo and predicted structures.
- Monomer input structures will be sampled from paired structures available in PINDER/PLINDER, balanced between apo and predicted structures and stratified by "flexibility" level according to specified conformational difference thresholds.

## Winner’s presentation

The winners will be invited to present their work at the MLSB workshop.
