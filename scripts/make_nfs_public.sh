#!/bin/bash

if [ "$#" -lt 2 ]; then
  echo "Usage: $0 <PLINDER_RELEASE> <PLINDER_ITERATION>"
  exit 1
fi

PLINDER_RELEASE="$1"
PLINDER_ITERATION="$2"

gsutil -m cp -r /plinder/${PLINDER_RELEASE}/index/* gs://plinder/${PLINDER_RELEASE}/${PLINDER_ITERATION}/index/
gsutil -m cp -r /plinder/${PLINDER_RELEASE}/clusters/* gs://plinder/${PLINDER_RELEASE}/${PLINDER_ITERATION}/clusters/
gsutil -m cp -r /plinder/${PLINDER_RELEASE}/entries/* gs://plinder/${PLINDER_RELEASE}/${PLINDER_ITERATION}/entries/
gsutil -m cp -r /plinder/${PLINDER_RELEASE}/fingerprints/* gs://plinder/${PLINDER_RELEASE}/${PLINDER_ITERATION}/fingerprints/
gsutil -m cp -r /plinder/${PLINDER_RELEASE}/ligand_scores/* gs://plinder/${PLINDER_RELEASE}/${PLINDER_ITERATION}/ligand_scores/
gsutil -m cp -r /plinder/${PLINDER_RELEASE}/ligands/* gs://plinder/${PLINDER_RELEASE}/${PLINDER_ITERATION}/ligands/
gsutil -m cp -r /plinder/${PLINDER_RELEASE}/mmp/* gs://plinder/${PLINDER_RELEASE}/${PLINDER_ITERATION}/mmp/
gsutil -m cp -r /plinder/${PLINDER_RELEASE}/scores/* gs://plinder/${PLINDER_RELEASE}/${PLINDER_ITERATION}/scores/
gsutil -m cp -r /plinder/${PLINDER_RELEASE}/archives/* gs://plinder/${PLINDER_RELEASE}/${PLINDER_ITERATION}/systems/
gsutil -m cp -r /plinder/${PLINDER_RELEASE}/splits/* gs://plinder/${PLINDER_RELEASE}/${PLINDER_ITERATION}/splits/
gsutil -m cp -r /plinder/${PLINDER_RELEASE}/leakage/* gs://plinder/${PLINDER_RELEASE}/${PLINDER_ITERATION}/leakage/
gsutil -m cp -r /plinder/${PLINDER_RELEASE}/linked_systems/* gs://plinder/${PLINDER_RELEASE}/${PLINDER_ITERATION}/linked_systems/
gsutil -m cp -r /plinder/${PLINDER_RELEASE}/dbs/subdbs/apo_foldseek/ids.tsv gs://plinder/${PLINDER_RELEASE}/${PLINDER_ITERATION}/dbs/subdbs/apo_foldseek/ids.tsv
gsutil -m cp -r /plinder/${PLINDER_RELEASE}/dbs/subdbs/holo_foldseek/ids.tsv gs://plinder/${PLINDER_RELEASE}/${PLINDER_ITERATION}/dbs/subdbs/holo_foldseek/ids.tsv
gsutil -m cp -r /plinder/${PLINDER_RELEASE}/dbs/subdbs/pred_foldseek/ids.tsv gs://plinder/${PLINDER_RELEASE}/${PLINDER_ITERATION}/dbs/subdbs/pred_foldseek/ids.tsv
gsutil -m cp -r /plinder/${PLINDER_RELEASE}/dbs/subdbs/apo_mmseqs/ids.tsv gs://plinder/${PLINDER_RELEASE}/${PLINDER_ITERATION}/dbs/subdbs/apo_mmseqs/ids.tsv
gsutil -m cp -r /plinder/${PLINDER_RELEASE}/dbs/subdbs/holo_mmseqs/ids.tsv gs://plinder/${PLINDER_RELEASE}/${PLINDER_ITERATION}/dbs/subdbs/holo_mmseqs/ids.tsv
gsutil -m cp -r /plinder/${PLINDER_RELEASE}/dbs/subdbs/pred_mmseqs/ids.tsv gs://plinder/${PLINDER_RELEASE}/${PLINDER_ITERATION}/dbs/subdbs/pred_mmseqs/ids.tsv
