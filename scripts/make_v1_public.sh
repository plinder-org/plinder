#!/bin/bash

# datasets
gsutil -m cp -r gs://plinder-collab-bucket/2024-04/v1/clusters/* gs://plinder/2024-04/v1/clusters/
gsutil -m cp -r gs://plinder-collab-bucket/2024-04/v1/entries/* gs://plinder/2024-04/v1/entries/
gsutil -m cp -r gs://plinder-collab-bucket/2024-04/v1/fingerprints/* gs://plinder/2024-04/v1/fingerprints/
gsutil -m cp -r gs://plinder-collab-bucket/2024-04/v1/ligand_scores/* gs://plinder/2024-04/v1/ligand_scores/
gsutil -m cp -r gs://plinder-collab-bucket/2024-04/v1/ligands/* gs://plinder/2024-04/v1/ligands/
gsutil -m cp -r gs://plinder-collab-bucket/2024-04/v1/mmp/* gs://plinder/2024-04/v1/mmp/
gsutil -m cp -r gs://plinder-collab-bucket/2024-04/v1/scores/* gs://plinder/2024-04/v1/scores/
gsutil -m cp -r gs://plinder-collab-bucket/2024-04/v1/systems/* gs://plinder/2024-04/v1/systems/
# don't bring index cruft
gsutil cp gs://plinder-collab-bucket/2024-04/v1/index/annotation_table.parquet gs://plinder/2024-04/v1/index/
gsutil cp gs://plinder-collab-bucket/2024-04/v1/index/annotation_table_nonredundant.parquet gs://plinder/2024-04/v1/index/
# publishable splits
gsutil -m cp -r gs://plinder-collab-bucket/2024-04/v1/splits/batch_6/* gs://plinder/2024-04/v1/splits/batch_6/
gsutil -m cp -r gs://plinder-collab-bucket/metaflow/config/split/batch_6/* gs://plinder/2024-04/v1/splits/batch_6/
gsutil cp gs://plinder-collab-bucket/2024-04/v1/splits/v0/splits/plinder_split_v0.parquet gs://plinder/2024-04/v1/splits/v0/plinder-v0.parquet
gsutil cp gs://plinder-collab-bucket/2024-04/v1/splits/v0/splits/plinder_split_v0.yaml gs://plinder/2024-04/v1/splits/v0/plinder-v0.yaml
gsutil cp gs://plinder-collab-bucket/2024-04/v1/workshop_split/ecod_split.csv gs://plinder/2024-04/v1/splits/v0/plinder-ecod.csv
gsutil cp gs://plinder-collab-bucket/2024-04/v1/workshop_split/time_split.csv gs://plinder/2024-04/v1/splits/v0/plinder-time.csv
gsutil cp gs://plinder-collab-bucket/2024-04/v1/workshop_split/plinder_v0_no_posbuster.csv gs://plinder/2024-04/v1/splits/v0/plinder-no-posebusters.csv
# only copy leakage for publishable splits
gsutil -m cp -r gs://plinder-collab-bucket/2024-04/v1/leakage/plinder* gs://plinder/2024-04/v1/leakage/
gsutil -m cp -r gs://plinder-collab-bucket/2024-04/v1/leakage/splits* gs://plinder/2024-04/v1/leakage/
