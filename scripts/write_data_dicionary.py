#!/usr/bin/env python
from pathlib import Path

import pandas as pd

# TODO: is this function used anywhere?
df = pd.concat([
    pd.read_csv(path, sep='\t') for path in
    Path('src/plinder-data/column_descriptions').rglob('*.tsv')
]).reset_index(drop=True)
with open('src/plinder-data/data_dictionary.md', 'w') as f:
    f.write(df.to_markdown() + '\n')
