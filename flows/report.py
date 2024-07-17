# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
import os
import json
from dataclasses import dataclass, asdict
from pathlib import Path


@dataclass
class PDB:
    two_char_code: str
    pdb_id: str
    ingest: int
    report: int
    entry: int
    num_holo: int
    num_apo: int
    num_system: int


def get_pdb_data(base: str, pdb_dir: str) -> PDB:
    two_char_code = pdb_dir[-3:-1]
    pdb_id = pdb_dir[-4:]
    entry_dir = f"{base}/entries/{two_char_code}"
    num_holo = 0
    num_apo = 0
    num_system = 0
    entry = int(os.path.exists(f"{entry_dir}/{pdb_id}.json"))
    if entry:
        with open(f"{entry_dir}/{pdb_id}.json") as f:
            contents = json.load(f)
            num_system = len(contents["systems"])
            for chain in contents.get("chains", {}).values():
                if chain.get("holo"):
                    num_holo += 1
                else:
                    num_apo += 1
    return PDB(
        two_char_code=two_char_code,
        pdb_id=pdb_id,
        ingest=1,
        report=int(os.path.exists(f"{base}/reports/{two_char_code}/{pdb_id}/{pdb_id}_validation.xml.gz")),
        entry=entry,
        num_holo=num_holo,
        num_apo=num_apo,
        num_system=num_system,
    )


def main(upload: bool = False):
    import pandas as pd

    plinder_mount = os.getenv("PLINDER_MOUNT")
    plinder_release = os.getenv("PLINDER_RELEASE")
    base = f"{plinder_mount}/{plinder_release}"
    total_cifs = 0
    total_reports = 0
    Path("reports").mkdir(exist_ok=True)
    pdbs = []
    fails = []
    for i, two_char_code in enumerate(sorted(os.listdir(f"{base}/ingest"))):
        if not i % 25:
            print(i, end=" ", flush=True)
        pdbs.extend([
            get_pdb_data(base, pdb_dir) for pdb_dir in
            os.listdir(f"{base}/ingest/{two_char_code}")
        ])
        fail_log = f"{base}/entries/{two_char_code}/failures.csv"
        if os.path.exists(fail_log):
            fails.append(pd.read_csv(fail_log))
    pdbs = pd.DataFrame([asdict(pdb) for pdb in pdbs])
    failures = pd.concat(fails).drop_duplicates(subset="pdb_id", keep="last").reset_index(drop=True)
    pdbs.to_csv("reports/pdbs.csv", index=False)
    failures.to_csv("reports/failures.csv", index=False)
    if upload:
        import gcs

        gcs.upload_report_from_file("reports/pdbs.csv")
        gcs.upload_report_from_file("reports/failures.csv")

    # print
    total_cifs = len(pdbs.index)
    total_reports = pdbs["report"].sum()
    total_entries = pdbs["entry"].sum()
    total_systems = pdbs["num_system"].sum()
    total_holo = pdbs["num_holo"].sum()
    total_apo = pdbs["num_apo"].sum()
    missing_entries = total_cifs - total_entries
    failed_entries = len(failures.index)
    msg_counts = failures["error"].value_counts()
    print(f"""\
total mmcifs          : {total_cifs}
total reports         : {total_reports}
total entries         : {total_entries}
total systems         : {total_systems}
total holo chains     : {total_holo}
total apo chains      : {total_apo}
missing entries       : {missing_entries}
failed entries        : {failed_entries}
most common failure   : {msg_counts.iloc[0]}

CSVs saved to reports/ directory for further analysis
""")
    return {
        "pdbs": pdbs,
        "failures": failures,
    }


if __name__ == "__main__":
    main()
