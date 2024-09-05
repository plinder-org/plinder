# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
from dataclasses import dataclass, field
from functools import partial
from hashlib import md5
from io import StringIO
from json import dumps
from os import getenv
from pathlib import Path
from typing import Any, Optional

from omegaconf import DictConfig, ListConfig, OmegaConf

from plinder.core.utils.log import setup_logger

LOG = setup_logger(__name__)


def _validate_cfg(*, cfg: DictConfig, schema: dict[str, Any]) -> DictConfig:
    """
    Reinitialize through dataclasses to get post-init
    validation logic from parameters that are only picked up through
    omegaconf (e.g. through the CLI).

    Parameters
    ----------
    cfg : DictConfig
        the omegaconf container
    schema : dict[str, Any]
        the schema to validate the config against

    Returns
    -------
    DictConfig
        the validated config with post-init validation logic
    """
    keys = set(cfg.keys()).union(set(schema.keys()))
    cfg = OmegaConf.to_container(cfg)
    cfg.get("data", {}).pop("plinder_dir", None)
    cfg.get("data", {}).pop("plinder_remote", None)
    return DictConfig({str(k): schema[str(k)](**cfg.get(k, {})) for k in keys})


def _clean_sort_config(*, cfg: Any) -> Any:
    """
    Recursively sort and clean a nested object to
    pure python types for consistent hashing.

    Parameters
    ----------
    cfg : Any
        object, omegaconf container, or primitive

    Returns
    -------
    Any
        the cleaned and sorted object (if nested)
    """
    if isinstance(cfg, (dict, DictConfig)):
        return {k: _clean_sort_config(cfg=v) for k, v in sorted(cfg.items())}
    elif isinstance(cfg, (list, ListConfig)):
        return [_clean_sort_config(cfg=v) for v in sorted(cfg)]
    else:
        return cfg


class _get_config:
    _schema: dict[str, Any] = {}
    _packages: dict[str, set[str]] = {}
    _cfg = DictConfig({})

    def _clear(self) -> None:
        self._schema = {}
        self._packages = {}
        self._cfg = DictConfig({})

    def __call__(
        self,
        *,
        schema: dict[str, Any],
        package_schema: str,
        config: Optional[dict[str, Any]] = None,
        config_file: Optional[str] = None,
        config_contents: Optional[str] = None,
        config_args: Optional[list[str]] = None,
        cached: bool = True,
    ) -> DictConfig:
        """
        Start with the default config, then add any
        config loaded from a yaml file (if provided)
        either as a path or a string and then pull in
        any CLI arguments passed.

        Parameters
        ----------
        config: dict[str, Any], default=None
            config overrides to use
        config_file : str, default=None
            path to a yaml file
        config_contents : str, default=None
            string contents of a yaml file
        config_args : list[str], default=None
            provided mainly to support testing CLI overrides
        schema : dict[str, Any], default=None
            schema to validate the config against
        package_schema: str
            the plinder subpackage name
        cached : bool, default=True
            if True, return cached config if it exists

        Returns
        -------
        config : DictConfig
            the fully resolved, merged config
        """
        if (
            cached
            and len(self._cfg)
            and schema.items() <= self._schema.items()
            and config is None
            and config_file is None
            and config_contents is None
            and config_args is None
        ):
            return self._cfg
        if package_schema not in self._packages:
            dups = set(self._schema.keys()).intersection(schema.keys())
            degenerate = False
            for keys in self._packages.values():
                if dups <= keys:
                    degenerate = True
                    break
            if len(dups) and not degenerate:
                raise ValueError(
                    f"schema keys must be unique: found keys {dups} in schema='{package_schema}'"
                )
            self._packages[package_schema] = set(schema.keys())
        self._schema.update(schema)
        args = []
        if config is not None:
            args.append(DictConfig(config))
        if config_file is not None:
            # load config from a yaml file
            args.append(DictConfig(OmegaConf.load(config_file)))
        if config_contents is not None:
            # load config from a yaml string
            LOG.debug("loading config from string, will skip CLI overrides")
            args.append(DictConfig(OmegaConf.load(StringIO(config_contents))))
        try:
            # catchall for non-standard execution models
            # e.g. in metaflow or jupyter
            cli = OmegaConf.from_cli(config_args)
            # omegaconf casts single int-like two_char_codes to int
            if cli.get("scatter", {}).get("two_char_codes") is not None:
                cli.scatter.two_char_codes = str(cli.scatter.two_char_codes)
            if len(cli):
                _validate_cfg(cfg=cli, schema=self._schema)
                args.append(cli)
        except Exception as e:
            LOG.debug(f"failed to validate cli args: {e}")
        cfg = DictConfig({})
        if len(args):
            cfg = DictConfig(OmegaConf.merge(*args))
        cfg = _validate_cfg(cfg=cfg, schema=self._schema)
        if cached:
            self._cfg = cfg
        return cfg


def get_config_hash(config_obj: Any) -> str:
    """
    Return a consistent hash of any object, dataclass, or omegaconf container.

    Parameters
    ----------
    config_obj : dict | dataclass | DictConfig
        the configuration to hash

    Returns
    -------
    config_hash : str
        the hash of the configuration
    """
    if isinstance(config_obj, (dict, DictConfig)):
        contents = config_obj
    else:
        contents = config_obj.__dict__
    return md5(dumps(_clean_sort_config(cfg=contents)).encode("utf-8")).hexdigest()


def _getenv_default(key: str, default: str) -> str:
    return getenv(key, default)


@dataclass
class DataConfig:
    """
    A class for all the data configuration.

    Notes
    -----
    plinder_dir is always derived from plinder_mount/plinder_release
    so do not set it manually!

    Attributes
    ----------
    plinder_release : str
        the plinder dataset version
    plinder_iteration : str
        the plinder dataset iteration
    plinder_mount : str, default="~/.local/share/plinder"
        the resting place for the plinder dataset
    plinder_bucket : str, default="plinder"
        the plinder bucket
    plinder_dir : str
        set automatically
    plinder_remote : str
        set automatically
    """

    plinder_release: str = field(
        default_factory=partial(_getenv_default, "PLINDER_RELEASE", "2024-06")
    )
    plinder_iteration: str = field(
        default_factory=partial(_getenv_default, "PLINDER_ITERATION", "v2")
    )
    plinder_mount: str = field(
        default_factory=partial(
            _getenv_default, "PLINDER_MOUNT", (Path.home() / ".local/share").as_posix()
        )
    )
    plinder_bucket: str = field(
        default_factory=partial(_getenv_default, "PLINDER_BUCKET", "plinder")
    )
    plinder_dir: str = field(init=False)
    plinder_remote: str = field(init=False)

    ingest: str = "ingest"
    validation: str = "validation"
    clusters: str = "clusters"
    entries: str = "entries"
    fingerprints: str = "fingerprints"
    fingerprint_file: str = "ligands_per_system.parquet"
    index: str = "index"
    ligand_scores: str = "ligand_scores"
    ligands: str = "ligands"
    links: str = "links"
    linked_structures: str = "linked_structures"
    mmp: str = "mmp"
    scores: str = "scores"
    splits: str = "splits"
    split_file: str = "split.parquet"
    systems: str = "systems"
    index_file: str = "annotation_table.parquet"
    force_update: bool = False

    def __post_init__(self) -> None:
        suffix = self.plinder_release
        if self.plinder_iteration:
            suffix = f"{self.plinder_release}/{self.plinder_iteration}"
        if self.plinder_mount in ["/plinder", "/", ""]:
            self.plinder_dir = f"{self.plinder_mount}/{suffix}"
        else:
            self.plinder_dir = f"{self.plinder_mount}/{self.plinder_bucket}/{suffix}"
        self.plinder_remote = f"gs://{self.plinder_bucket}/{suffix}"


@dataclass
class ContextConfig:
    two_char_codes: Any = ""
    pdb_ids: Any = ""
    system_ids: Any = ""

    def __post_init__(self) -> None:
        for attr in [
            "two_char_codes",
            "pdb_ids",
            "system_ids",
        ]:
            value = getattr(self, attr, None)
            if isinstance(value, str):
                setattr(self, attr, [val for val in value.split(",") if val])
            elif value is None:
                setattr(self, attr, [])


SCHEMA = {
    "data": DataConfig,
    "context": ContextConfig,
}

_config = _get_config()
get_config = partial(_config, schema=SCHEMA, package_schema="core")
