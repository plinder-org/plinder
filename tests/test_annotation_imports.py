"""Optional annotation dependencies are required only by their operations."""

import subprocess
import sys


def test_annotations_import_without_openstructure_or_networkit():
    subprocess.run(
        [
            sys.executable,
            "-c",
            """
import importlib.abc
import sys

class NoOptionalPackages(importlib.abc.MetaPathFinder):
    def find_spec(self, fullname, path=None, target=None):
        if fullname.split('.')[0] in {'ost', 'PDBValidation', 'networkit'}:
            raise ModuleNotFoundError(fullname)

sys.meta_path.insert(0, NoOptionalPackages())
from plinder.data.annotations.aggregate_annotations import Entry
from plinder.data.annotations.ligand_utils import Ligand
from plinder.data.annotations.protein_utils import Chain

for model in (Entry, Ligand, Chain):
    assert model.model_json_schema()['properties']
assert not {'ost', 'PDBValidation', 'networkit'} & sys.modules.keys()
""",
        ],
        check=True,
        capture_output=True,
        text=True,
    )
