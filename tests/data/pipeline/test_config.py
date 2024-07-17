# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
import sys
import unittest.mock
from textwrap import dedent
from omegaconf import OmegaConf, DictConfig

from plinder.data.pipeline import config

import pytest


@pytest.mark.parametrize("value, raises", [
    (0, True),
    (1, False),
    (2, False),
])
def test_foldseek_config(value, raises):
    if raises:
        with pytest.raises(ValueError):
            config.FoldseekConfig(alignment_type=value)
    else:
        config.FoldseekConfig(alignment_type=value)


def test_ingest_config():
    dc = config.IngestConfig()
    cfg = OmegaConf.structured(config.IngestConfig())
    assert dc.plinder_mount == cfg.plinder_mount


def test_default_config():
    cfg = config.get_config(cached=False)
    assert cfg.ingest.plinder_release is not None


def test_get_config_metaflow(tmp_path):
    file = tmp_path / "conf.yaml"
    file.write_text(
        dedent(
            """
            ingest:
              skip_specific_stages: foo
            """
        )
    )
    contents = dedent(
        """
        scatter:
          two_char_codes: xx
        """
    )
    cfg = config.get_config(
        cached=False,
        config_file=file.as_posix(),
        config_contents=contents,
    )
    assert cfg.ingest.skip_specific_stages == ["foo"]
    assert cfg.scatter.two_char_codes == ["xx"]


def test_get_config_comma_delimited():
    contents = dedent(
        """
        scatter:
          two_char_codes: xx,yy,zz
        """
    )
    cfg = config.get_config(
        cached=False,
        config_contents=contents,
    )
    assert cfg.scatter.two_char_codes == ["xx", "yy", "zz"]


def test_get_config_cli():
    test_args = ["prog", "scatter.two_char_batch_size=4"]
    with unittest.mock.patch("sys.argv", test_args):
        cfg = config.get_config(cached=False)
        assert cfg.scatter.two_char_batch_size == 4
