# Changelog

- WIP (Current — unreleased):
    - **Public downloads from Cameo**: releases download directly from `https://cameo3d.org/plinder/PLINDER-<YYYY-MM>/` (currently `2026-09`); `PLINDER_MIRROR_URL` selects another HTTP file server with the same directory layout. Downloads stream to disk, validate HTTP Content-Length, and replace cached files atomically. Interrupted transfers retry from the start. Cached files are reused when their size and Last-Modified time match the server, so data hotfixes are fetched on the next access. Google storage dependencies are optional: install `plinder[data]` for the GCS data-generation utilities.
    - **Major backend refactor**: replaced OST, gemmi, plip, openbabel with biotite + peppr for data generation; removed 6 dependencies from ingest pipeline
    - **Nucleic acid support**: DNA/RNA chains now correctly included as receptor neighbors, mainchain/sidechain detection works for both protein and nucleic acids ([#61](https://github.com/plinder-org/plinder/issues/61))
    - **Custom CIF support**: new `Entry.from_custom_cif_file` for structure-prediction outputs (Boltz, AlphaFold3, Chai-1) that ship CIFs without `_chem_comp_bond` ([#117](https://github.com/plinder-org/plinder/issues/117)). Bond orders come from `ligand_smiles_dict` via positional atom-order match (the convention these tools follow); element/count mismatches raise with the offending position, `force_substructure_match=True` opts into substructure matching when atom order isn't preserved. User SMILES win over CCD for both `smiles` and `resolved_stereo_matches_template` — closes a silent gap where biotite's `LIG` placeholder would pass any 3D conformer. Input CIFs are never mutated; optional `save_fixed_cif` persists the enriched copy.
    - **Stricter CIF ingest**: H/D filtered consistently (biotite's `filter_heavy`); multi-model CIFs warn and use model 1; multi-instance custom comp_ids must share heavy-atom naming (since `_chem_comp_bond` is comp_id-keyed); silent `connect_via_residue_names` and half-sanitized substructure fallbacks replaced with `ValueError` so corrupt inputs fail loudly.
    - **Stereochemistry**: CCD ideal 3D coordinates used as stereo ground truth; new `resolved_stereo_matches_template` flag validates resolved structure chirality against CCD template (handles partial resolution via MCS trimming)
    - **Interactions**: water bridge and metal bridge detection via peppr; halogen bond sidechain flag now computed (was hardcoded)
    - **Binding affinity**: fixed BindingDB matching — target sequence now validated against PDB SEQRES with 100% core identity, terminal tags/truncations tolerated ([#94](https://github.com/plinder-org/plinder/issues/94)); updated code to get the latest BindingDB release
    - **Optional eval**: `pip install plinder[eval]` adds PoseBusters ligand validation; OpenStructure-backed metrics use the Conda-only `openstructure` package; PoseBusters no longer runs during ingest
    - **PlinderSystem API**: new `receptor_structure` and `ligand_structures` (Biotite AtomArray) plus `ligand_mols` (RDKit Mol) properties; OpenStructure is confined to the evaluation implementation
    - **Chain type support**: `Chain.from_cif_data` now assigns proper one-letter codes and chem_types for nucleotides (`RNA Linking`, `DNA Linking`); new `Residue.is_modified` property covers both protein PTMs and modified nucleotide bases
    - **Save utils**: receptor/ligand chain naming generalized (`PDB_RECEPTOR_CHAINS`); system saving works for protein, NA, and mixed complexes
    - **System definition**: unified `min_polymer_size=12` replaces separate `min_polymer_size`/`max_non_small_mol_ligand_length` — polymers ≥ 12 residues are receptor, shorter are ligands (threshold matches minimum MMseqs2/Foldseek search length); molecules with BIRD annotation are ligands irrespective of size; ligand chains no longer appear in both receptor and ligand parts of system IDs.
    - **System grouping**: pocket-based grouping (≥ 3 shared receptor residues on the same chain instance) for adjacent binding sites (e.g. orthosteric + allosteric, cofactor + substrate in same active site); artifacts attach only via 4 Å proximity
    - **Dead code removal**: removed unused OST-based functions, PDB string roundtrips, duplicate SMILES derivation paths, v1 template matching (consolidated to Rascal MCES `get_matched_template`)
    - **License**: changed from GPL-2.0 to Apache-2.0 (GPL was only required by PLIP, now removed)

- 2024-06/v2:
    - New systems added based on the 2024-06 RCSB sync
    - Updated system definition to be more stable and depend only on ligand distance rather than PLIP
    - Added annotations for crystal contacts
    - Improved ligand handling and saving to fix some bond order issues
    - Improved covalency detection and annotation to reference each bond explicitly
    - Added linked apo/pred structures to v2/links and v2/linked_structures
    - <del>Added binding affinity annotations from [BindingDB](https://bindingdb.org)</del> (see known bugs!)
    - Added statistics requirement and other changes in the split to enrich test set diversity

- 2024-04/v1: Version described in the preprint, with updated redundancy removal by protein pocket and ligand similarity.
- 2024-04/v0: Version used to re-train DiffDock in the paper, with redundancy removal based on \<pdbid\>\_\<ligand ccd codes\>
