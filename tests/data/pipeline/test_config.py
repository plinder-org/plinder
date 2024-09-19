# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
import unittest.mock
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


def test_flow_config():
    dc = config._config.DataConfig()
    cfg = OmegaConf.structured(config._config.DataConfig())
    assert dc.plinder_mount == cfg.plinder_mount


def test_default_config():
    cfg = config.get_config(cached=False)
    assert cfg.data.plinder_release is not None


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
