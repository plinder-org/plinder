# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
from dataclasses import asdict
from hashlib import md5
from json import dumps

from plinder.core.utils import config


def test_get_config():
    cfg = config.get_config(config_args=[], cached=False)
    ocfg = config.DataConfig()
    ccfg = config.ContextConfig()
    assert cfg == {"data": asdict(ocfg), "context": asdict(ccfg)}


def test_get_config_passed():
    cfg = config.get_config(
        config_args=[], config={"data": {"plinder_release": "test"}}, cached=False
    )
    assert cfg.data.plinder_release == "test"


def test_get_config_env(monkeypatch):
    # dataclasses defer getenv to init so should be in sync
    monkeypatch.setenv("PLINDER_RELEASE", "test")
    cfg = config.get_config(config_args=[], cached=False)
    ocfg = config.DataConfig()
    assert cfg == {"data": asdict(ocfg), "context": asdict(config.ContextConfig())}
    assert ocfg.plinder_release == "test"
    assert cfg.data.plinder_release == "test"


def test_get_config_cli():
    # dataclasses don't support cli overrides
    cfg = config.get_config(config_args=["data.plinder_release=test"], cached=False)
    assert cfg.data.plinder_release == "test"


def test_get_config_contents():
    contents = """
data:
    plinder_release: test
"""
    cfg = config.get_config(config_contents=contents, config_args=[], cached=False)
    assert cfg.data.plinder_release == "test"


def test_get_config_file(tmp_path):
    file = tmp_path / "test.yaml"
    file.write_text(
        """
data:
    plinder_release: test
"""
    )
    cfg = config.get_config(
        config_file=str(file.as_posix()), config_args=[], cached=False
    )
    assert cfg.data.plinder_release == "test"


def test_get_config_hash():
    cfg = config.get_config(config_args=[], cached=False)
    assert config.get_config_hash(cfg)
    dc = config.DataConfig()
    assert config.get_config_hash(cfg.data) == config.get_config_hash(dc)
    dct = asdict(dc)
    assert config.get_config_hash(cfg.data) == config.get_config_hash(dct)


def test_get_config_hash_sorted():
    a = {"a": 1, "b": 2}
    b = {"b": 2, "a": 1}
    assert config.get_config_hash(a) == config.get_config_hash(b)
    assert (
        md5(dumps(a).encode("utf-8")).hexdigest()
        != md5(dumps(b).encode("utf-8")).hexdigest()
    )


def test_get_config_hash_sorted_nested():
    a = {"a": 1, "b": 2, "c": {"d": 3, "e": 4}}
    b = {"b": 2, "a": 1, "c": {"e": 4, "d": 3}}
    assert config.get_config_hash(a) == config.get_config_hash(b)
    assert (
        md5(dumps(a).encode("utf-8")).hexdigest()
        != md5(dumps(b).encode("utf-8")).hexdigest()
    )


def test_get_config_hash_sorted_nested_list():
    a = {"a": 1, "b": 2, "c": {"d": 3, "e": 4, "f": [1, 2, 3]}}
    b = {"b": 2, "a": 1, "c": {"e": 4, "d": 3, "f": [3, 2, 1]}}
    assert config.get_config_hash(a) == config.get_config_hash(b)
    assert (
        md5(dumps(a).encode("utf-8")).hexdigest()
        != md5(dumps(b).encode("utf-8")).hexdigest()
    )
