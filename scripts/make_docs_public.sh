#!/bin/bash

cp src/plinder-data/README.md src/plinder-data/README_v2.md
sed -i.tmp 's/2024-04/2024-06/g' src/plinder-data/README_v2.md
sed -i.tmp 's/v1/v2/g' src/plinder-data/README_v2.md
rm src/plinder-data/README_v2.md.tmp

echo "remember to upload them! see the rest of this script"
echo "add the following to README_v2.md

    |-- linked_structures
    |   |-- {search_db}_links.parquet


|   |-- linked_structures    # linked structures
"

python -c "from pathlib import Path; import pandas as pd; df = pd.concat([pd.read_csv(path, sep='\t') for path in Path('src/plinder-data/column_descriptions').rglob('*.tsv')]).reset_index(drop=True); f = open('src/plinder-data/data_dictionary.md', 'w'); f.write(df.to_markdown() + '\n'); f.close()"
#
# gsutil cp src/plinder-data/data_dictionary.md gs://plinder/2024-04/v1/data_dictionary.md
# gsutil cp src/plinder-data/manifest.md gs://plinder/manifest.md
# gsutil cp src/plinder-data/README.md gs://plinder/2024-04/v1/README.md
# gsutil cp src/plinder-data/README_v2.md gs://plinder/2024-04/v1/README.md
