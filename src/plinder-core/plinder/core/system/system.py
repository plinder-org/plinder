from typing import Any, Optional

from omegaconf import DictConfig

from plinder.core.index import utils
from plinder.core.utils.log import setup_logger

LOG = setup_logger(__name__)


class PlinderSystem:
    def __init__(
        self,
        system_id: str,
        entry: dict[str, Any],
    ) -> None:
        self.system_id = system_id
        self._entry = entry
        if system_id not in entry["systems"]:
            raise ValueError(f"system_id={system_id} not in entry")
        self.system = entry["systems"][system_id]

    def structures(self) -> list[str]:
        return []

    @classmethod
    def from_system_id(
        cls,
        system_id: str,
        cfg: Optional[DictConfig] = None,
        prune: bool = True,
    ) -> Optional["PlinderSystem"]:
        man = utils.get_manifest(cfg=cfg)
        if system_id not in man["system_id"].values:
            LOG.error(f"system_id={system_id} not in manifest")
            return None
        entry_pdb_id = system_id.split("__")[0]
        entry = utils.load_entries(pdb_ids=[entry_pdb_id], cfg=cfg, prune=prune)
        system = entry[entry_pdb_id]["systems"].get(system_id)
        if system is None:
            LOG.error(f"system_id={system_id} not in entry={entry_pdb_id}")
            return None
        return cls(
            system_id=system_id,
            entry=entry,
        )
