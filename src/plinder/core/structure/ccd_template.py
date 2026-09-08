# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
"""CCD component templates: missing-atom detection and opt-in completion."""
from __future__ import annotations

from collections.abc import Iterable
from functools import cache
from typing import NamedTuple

import biotite.structure as struc
import biotite.structure.info as bt_info
import numpy as np
from biotite.structure.atoms import AtomArray

from plinder.core.utils.log import setup_logger

LOG = setup_logger(__name__)

UNRESOLVED_ANNOTATION = "is_unresolved"


class ComponentTemplate(NamedTuple):
    heavy: tuple[str, ...]  # heavy atom names in CCD order
    elements: dict[str, str]
    leaving: frozenset[str]
    bonds: tuple[tuple[str, str, int], ...]  # heavy-heavy bonds by atom name


@cache
def ccd_component_template(comp_id: str) -> ComponentTemplate | None:
    """Heavy atoms, elements, leaving atoms and bonds of a CCD component."""
    from biotite.structure.info import ccd as bundled_ccd

    try:
        template = bt_info.residue(comp_id, allow_missing_coord=True)
        chem_comp_atom = bundled_ccd.get_ccd()["chem_comp_atom"]
    except Exception as exc:
        LOG.warning(f"CCD lookup failed for {comp_id}: {exc}")
        return None
    if template is None:
        return None
    heavy_mask = ~np.isin(template.element, ["H", "D"])
    heavy = tuple(str(name) for name in template.atom_name[heavy_mask])
    elements = {
        str(name): str(element)
        for name, element in zip(template.atom_name, template.element)
    }
    rows = chem_comp_atom["comp_id"].as_array(str) == comp_id
    leaving = frozenset(
        chem_comp_atom["atom_id"].as_array(str)[rows][
            chem_comp_atom["pdbx_leaving_atom_flag"].as_array(str)[rows] == "Y"
        ]
    )
    bonds: list[tuple[str, str, int]] = []
    if template.bonds is not None:
        heavy_set = set(heavy)
        for index_a, index_b, order in template.bonds.as_array():
            name_a, name_b = (
                str(template.atom_name[index_a]),
                str(template.atom_name[index_b]),
            )
            if name_a in heavy_set and name_b in heavy_set:
                bonds.append((name_a, name_b, int(order)))
    return ComponentTemplate(heavy, elements, leaving, tuple(bonds))


def ccd_heavy_atom_names(comp_id: str) -> tuple[frozenset[str], frozenset[str]] | None:
    """Heavy atom names and leaving-group atom names of a CCD component."""
    template = ccd_component_template(comp_id)
    if template is None:
        return None
    return frozenset(template.heavy), template.leaving


def unresolved_atoms_from_template(
    comp_id: str, resolved_atom_names: Iterable[str]
) -> list[str] | None:
    """CCD heavy atoms (leaving groups excluded) absent from ``resolved_atom_names``.

    Returns ``None`` when ``comp_id`` is not in the CCD.

    TODO: leaving atoms are excluded unconditionally; for a residue that is not
    covalently linked (free reducing-end sugar, C-terminal OXT) they are real
    atoms and their absence should count. Needs the struct_conn links per residue.
    """
    template = ccd_component_template(comp_id)
    if template is None:
        return None
    resolved = set(resolved_atom_names)
    return [
        name
        for name in template.heavy
        if name not in template.leaving and name not in resolved
    ]


def add_missing_atoms(
    atoms: AtomArray, *, mask_annotation: str = UNRESOLVED_ANNOTATION
) -> AtomArray:
    """Append CCD heavy atoms missing from each residue, with NaN coordinates.

    A residue is completed only when its CCD template is known and every
    resolved heavy atom name is a template name; leaving atoms and hydrogens are
    never added and fully unresolved residues are not created. Added atoms copy
    their residue's annotations, take ``element`` from the CCD, get
    ``occupancy``/``b_factor``/``charge`` 0 where present, and ``mask_annotation``
    True (deposited atoms get False). Existing bonds are kept and the template
    bonds of added atoms appended. The input is not modified.
    """
    count = atoms.array_length()
    starts = struc.get_residue_starts(atoms, add_exclusive_stop=True)
    solvent = struc.filter_solvent(atoms)
    hydrogen = np.isin(atoms.element, ["H", "D"])
    completed: list[tuple[int, int, list[str], ComponentTemplate]] = []
    mismatched = 0
    for start, stop in zip(starts[:-1], starts[1:]):
        if solvent[start]:
            continue
        template = ccd_component_template(str(atoms.res_name[start]))
        if template is None:
            continue
        resolved = {
            str(name)
            for name, is_h in zip(atoms.atom_name[start:stop], hydrogen[start:stop])
            if not is_h
        }
        if not resolved <= set(template.heavy):
            mismatched += 1
            continue
        missing = [
            name
            for name in template.heavy
            if name not in template.leaving and name not in resolved
        ]
        if missing:
            completed.append((int(start), int(stop), missing, template))

    categories = atoms.get_annotation_categories()
    existing_mask = (
        atoms.get_annotation(mask_annotation).astype(bool)
        if mask_annotation in categories
        else np.zeros(count, dtype=bool)
    )
    if not completed:
        result = atoms.copy()
        result.set_annotation(mask_annotation, existing_mask)
        if mismatched:
            LOG.info(
                "add_missing_atoms: %d residues skipped, names off-template", mismatched
            )
        return result

    positions = np.concatenate(
        [np.full(len(missing), stop) for _, stop, missing, _ in completed]
    )
    sources = np.concatenate(
        [np.full(len(missing), start) for start, _, missing, _ in completed]
    )
    new_names = np.array(
        [name for _, _, missing, _ in completed for name in missing], dtype=str
    )
    new_elements = np.array(
        [
            template.elements[name]
            for _, _, missing, template in completed
            for name in missing
        ],
        dtype=str,
    )
    added = len(positions)
    result = struc.AtomArray(count + added)
    result.coord = np.insert(
        atoms.coord,
        positions,
        np.full((added, 3), np.nan, dtype=atoms.coord.dtype),
        axis=0,
    )
    for category in categories:
        values = atoms.get_annotation(category)
        if category == "atom_name":
            fill = new_names
        elif category == "element":
            fill = new_elements
        elif category in {"occupancy", "b_factor", "charge"}:
            fill = np.zeros(added, dtype=values.dtype)
        elif category == mask_annotation:
            continue
        else:
            fill = values[sources]
        result.set_annotation(category, np.insert(values, positions, fill))
    result.set_annotation(
        mask_annotation, np.insert(existing_mask, positions, np.ones(added, dtype=bool))
    )

    if atoms.bonds is not None:
        new_index = np.arange(count) + np.searchsorted(
            positions, np.arange(count), "right"
        )
        bonds = atoms.bonds.as_array()
        remapped = bonds.copy()
        remapped[:, 0] = new_index[bonds[:, 0]]
        remapped[:, 1] = new_index[bonds[:, 1]]
        extra: list[list[int]] = []
        for start, stop, missing, template in completed:
            index_by_name = {
                str(name): int(new_index[i])
                for i, name in zip(range(start, stop), atoms.atom_name[start:stop])
            }
            first_new = int(new_index[stop - 1]) + 1
            index_by_name.update(
                {name: first_new + k for k, name in enumerate(missing)}
            )
            for name_a, name_b, order in template.bonds:
                if (
                    (name_a in missing or name_b in missing)
                    and name_a in index_by_name
                    and name_b in index_by_name
                ):
                    extra.append([index_by_name[name_a], index_by_name[name_b], order])
        all_bonds = (
            np.vstack([remapped, np.array(extra, dtype=int)]) if extra else remapped
        )
        result.bonds = struc.BondList(count + added, all_bonds)

    LOG.info(
        "add_missing_atoms: added %d atoms to %d residues; %d residues skipped, names off-template",
        added,
        len(completed),
        mismatched,
    )
    return result
