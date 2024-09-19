# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
import pytest
from plinder.data.pipeline import transform


@pytest.mark.parametrize(
    "pdb_range",
    [
        "SM:9-163",
        "SM:9--163",
        "SM:-9--163",
        "SM:9-163L",
        "B:-3-715",
        "A:-4-131",
        "B:-1-75",
        "B:-4-1927",
    ],
)
def test_parse_pdb_range(pdb_range):
    assert len(transform.parse_pdb_range(pdb_range)) == 3


def test_transform_ecod_data(tmp_path):
    ecod_path = tmp_path / "ecod_raw.tsv"
    sample = """\
#/data/ecod/database_versions/v291/ecod.develop291.domains.txt
#ECOD version develop291
#Domain list version 1.6
#Grishin lab (http://prodata.swmed.edu/ecod)
#uid	ecod_domain_id	manual_rep	t_id	pdb	chain	pdb_range	seqid_range	unp_acc	arch_name	x_name	h_name	t_name	f_name	asm_status	ligand
000000267	e1udzA1	MANUAL_REP	1.1.1	1udz	A	A:203-381	A:4-182	P56690	beta barrels	"cradle loop barrel"	"RIFT-related"	"acid protease"	F_UNCLASSIFIED	NOT_DOMAIN_ASSEMBLY	NO_LIGANDS_4A
000023408	e1ileA4	AUTO_NONREP	1.1.1	1ile	A	A:203-381	A:203-381	P56690	beta barrels	"cradle loop barrel"	"RIFT-related"	"acid protease"	F_UNCLASSIFIED	NOT_DOMAIN_ASSEMBLY	NO_LIGANDS_4A
000023411	e1ue0B1	AUTO_NONREP	1.1.1	1ue0	B	B:203-381	B:4-182	P56690	beta barrels	"cradle loop barrel"	"RIFT-related"	"acid protease"	F_UNCLASSIFIED	NOT_DOMAIN_ASSEMBLY	NO_LIGANDS_4A
000158260	e1wk8A1	AUTO_NONREP	1.1.1	1wk8	A	A:201-382	A:7-188	P56690	beta barrels	"cradle loop barrel"	"RIFT-related"	"acid protease"	F_UNCLASSIFIED	NOT_DOMAIN_ASSEMBLY	NO_LIGANDS_4A
001842922	e5fofA1	AUTO_NONREP	1.1.1	5fof	A	A:258-340,A:363-567	A:29-103,A:126-330	B3L7I1	beta barrels	"cradle loop barrel"	"RIFT-related"	"acid protease"	F_UNCLASSIFIED	NOT_DOMAIN_ASSEMBLY	NO_LIGANDS_4A
001842923	e5fofB1	AUTO_NONREP	1.1.1	5fof	B	B:258-340,B:363-567	B:29-103,B:126-330	B3L7I1	beta barrels	"cradle loop barrel"	"RIFT-related"	"acid protease"	F_UNCLASSIFIED	NOT_DOMAIN_ASSEMBLY	NO_LIGANDS_4A
"""
    ecod_path.write_text(sample)
    df = transform.transform_ecod_data(raw_ecod_path=ecod_path)
    assert len(df.index) == 8
