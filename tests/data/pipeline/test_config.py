# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
import unittest.mock
from dataclasses import asdict
from textwrap import dedent

import pytest
from omegaconf import OmegaConf
from plinder.data.pipeline import config


@pytest.mark.parametrize(
    "value, raises",
    [
        (0, True),
        (1, False),
        (2, False),
    ],
)
def test_foldseek_config(value, raises):
    if raises:
        with pytest.raises(ValueError):
            config.FoldseekConfig(alignment_type=value)
    else:
        config.FoldseekConfig(alignment_type=value)


def test_search_defaults_retain_ten_thousand_hits() -> None:
    assert config.FoldseekConfig().max_seqs == 10_000
    assert config.FoldseekConfig().min_seq_id == 0.2
    assert config.MMSeqsConfig().max_seqs == 10_000
    assert config.MMSeqsConfig().min_seq_id == 0.2


@pytest.mark.parametrize("config_type", [config.FoldseekConfig, config.MMSeqsConfig])
def test_search_hit_limit_must_be_positive(config_type) -> None:
    with pytest.raises(ValueError, match="max_seqs must be positive"):
        config_type(max_seqs=0)


@pytest.mark.parametrize("config_type", [config.FoldseekConfig, config.MMSeqsConfig])
@pytest.mark.parametrize("min_seq_id", [-0.1, 1.1])
def test_search_minimum_sequence_identity_must_be_a_fraction(
    config_type, min_seq_id
) -> None:
    with pytest.raises(ValueError, match="min_seq_id must be in"):
        config_type(min_seq_id=min_seq_id)


def test_flow_config():
    dc = config._config.DataConfig()
    cfg = OmegaConf.structured(config._config.DataConfig())
    assert dc.plinder_mount == cfg.plinder_mount


def test_default_config():
    cfg = config.get_config(cached=False)
    assert cfg.data.plinder_release is not None
    assert cfg.scorer.minimum_threshold == 0.3
    assert cfg.scorer.max_alignment_rows_per_query == 5_000_000
    assert cfg.scorer.max_query_protein_chains == 30
    assert cfg.scorer.max_query_proper_ligand_chains == 30
    assert list(cfg.flow.cluster_thresholds) == [30, 50, 70, 90, 100]
    assert cfg.flow.component_reduction_metric_workers == 4
    assert cfg.interface.contact_radius == 10.0
    assert cfg.interface.min_chain_length == 12
    assert cfg.interface.min_interface_residues == 7


def test_ingest_annotation_defaults_have_one_owner():
    assert asdict(config.AnnotationConfig()) == {
        "neighboring_residue_threshold": 6.0,
        "neighboring_ligand_threshold": 4.0,
        "min_polymer_size": 12,
        "min_shared_pocket_members": 3,
    }
    assert asdict(config.EntryConfig()) == {
        "interaction_search_threshold": 10.0,
        "tessellation_atom_limit": 2_000_000,
        "data_dir": None,
        "save_folder": None,
    }


@pytest.mark.parametrize("source", ["yaml", "cli"])
@pytest.mark.parametrize(
    ("include_ligands", "include_interfaces"),
    [(True, True), (True, False), (False, True)],
)
def test_annotation_overrides_reach_entry_reader(
    tmp_path, source, include_ligands, include_interfaces
):
    from plinder.data.get_system_annotations import GetPlinderAnnotation

    thresholds = {
        "neighboring_residue_threshold": 7.0,
        "neighboring_ligand_threshold": 5.0,
        "min_polymer_size": 15,
        "min_shared_pocket_members": 4,
    }
    if source == "yaml":
        contents = "annotation:\n" + "".join(
            f"  {name}: {value}\n" for name, value in thresholds.items()
        )
        contents += "entry:\n  interaction_search_threshold: 11.0\n"
        cfg = config.get_config(config_contents=contents, config_args=[], cached=False)
    else:
        args = [f"annotation.{name}={value}" for name, value in thresholds.items()]
        args.append("entry.interaction_search_threshold=11.0")
        cfg = config.get_config(config_args=args, cached=False)

    annotator = GetPlinderAnnotation(
        tmp_path / "source.cif",
        tmp_path / "validation.xml.gz",
        entry_cfg=dict(cfg.entry),
        **dict(cfg.annotation),
    )
    options = annotator._entry_options(
        include_ligands=include_ligands, include_interfaces=include_interfaces
    )
    assert {name: options[name] for name in thresholds} == thresholds
    assert options["interaction_search_threshold"] == 11.0
    assert options["include_ligands"] is include_ligands
    assert options["include_interfaces"] is include_interfaces


@pytest.mark.parametrize("name", list(asdict(config.AnnotationConfig())))
@pytest.mark.parametrize("source", ["yaml", "cli", "argv"])
def test_shared_thresholds_are_rejected_in_entry_section(name, source, monkeypatch):
    with pytest.raises(TypeError, match=name):
        if source == "yaml":
            config.get_config(
                config_contents=f"entry:\n  {name}: 5\n", config_args=[], cached=False
            )
        elif source == "cli":
            config.get_config(config_args=[f"entry.{name}=5"], cached=False)
        else:
            monkeypatch.setattr("sys.argv", ["ingest", f"entry.{name}=5"])
            config.get_config(cached=False)


@pytest.mark.parametrize(
    "kwargs",
    [
        {"contact_radius": 0},
        {"min_chain_length": 0},
        {"min_interface_residues": 0},
    ],
)
def test_interface_config_limits_must_be_positive(kwargs) -> None:
    with pytest.raises(ValueError, match="interface"):
        config.InterfaceConfig(**kwargs)


def test_component_reduction_metric_workers_must_be_positive() -> None:
    with pytest.raises(ValueError, match="component_reduction_metric_workers"):
        config.FlowConfig(component_reduction_metric_workers=0)


@pytest.mark.parametrize("metric", ["shape", "color", "sucos_shape"])
def test_raw_3d_ligand_diagnostics_cannot_be_clustered(metric) -> None:
    with pytest.raises(ValueError, match="sucos_shape_pocket_qcov"):
        config.get_config(
            cached=False,
            config={"flow": {"cluster_metrics": [metric]}},
        )


def test_alignment_mapping_row_budget_must_be_positive() -> None:
    with pytest.raises(ValueError, match="max_alignment_rows_per_query"):
        config.ScorerConfig(max_alignment_rows_per_query=0)


def test_get_config_metaflow(tmp_path):
    file = tmp_path / "conf.yaml"
    file.write_text(
        dedent(
            """
            flow:
              skip_specific_stages: foo
            """
        )
    )
    contents = dedent(
        """
        context:
          two_char_codes: xx
        """
    )
    cfg = config.get_config(
        cached=False,
        config_file=file.as_posix(),
        config_contents=contents,
    )
    assert cfg.flow.skip_specific_stages == ["foo"]
    assert cfg.context.two_char_codes == ["xx"]


def test_get_config_comma_delimited():
    contents = dedent(
        """
        context:
          two_char_codes: xx,yy,zz
        """
    )
    cfg = config.get_config(
        cached=False,
        config_contents=contents,
    )
    assert cfg.context.two_char_codes == ["xx", "yy", "zz"]


def test_get_config_cli():
    test_args = ["prog", "flow.download_rcsb_files_batch_size=4"]
    with unittest.mock.patch("sys.argv", test_args):
        cfg = config.get_config(cached=False)
        assert cfg.flow.download_rcsb_files_batch_size == 4
