# mypy: disable-error-code="attr-defined"

"""Vendored copy of peppr's ``sanitize`` with local over-valence fixes.

.. warning::

    TEMPORARY. This is a local, hand-maintained copy of the modified
    ``sanitize`` from unreleased peppr version. It carries boron-cage / main-group
    over-valence tolerance and related valence/kekulization fixes that are
    **not yet in a released** ``peppr`` and are not pip-installable from the
    internal repo.

    TODO(peppr): once these changes land in a published ``peppr`` release,
    DELETE this module and revert every import of the form
    ``from plinder.core.utils.sanitize import sanitize as peppr_sanitize``
    back to ``from peppr import sanitize as peppr_sanitize``
    (grep the tree for ``peppr_sanitize`` to find all call sites).

    Deliberate local deltas from the upstream peppr body: this warning docstring; the
    ``bool(...)`` wraps on the predicate returns (plinder's minimal mypy env has no rdkit
    stubs, so ``Chem`` is ``Any`` and ``strict``'s ``warn_return_any`` would flag them);
    and a trimmed ``_ORGANOMETALLIC_OMITTED_ELEMENTS`` comment (correctness, pending upstream).
"""

__all__ = ["sanitize"]

from collections.abc import Callable
from itertools import combinations, product

import rdkit
import rdkit.Chem.AllChem as Chem
from rdkit.Chem.rdmolops import SanitizeFlags

# Main-group metals/metalloids RDKit's SANITIZE_CLEANUP_ORGANOMETALLICS step leaves alone
# (it dativises only transition metals). An over-valent centre of one is tolerated neutral
# (:func:`_is_tolerated_over_valence`) and one counts as a coordinating metal when kekulizing
# (:func:`_is_metal`). Built from symbols for legibility, stored/compared as atomic numbers.
_ORGANOMETALLIC_OMITTED_ELEMENTS = frozenset(
    Chem.GetPeriodicTable().GetAtomicNumber(symbol)
    for symbol in (
        "B",
        "Li",
        "Be",
        "Na",
        "Mg",
        "Al",
        "K",
        "Ca",
        "Ga",
        "Ge",
        "As",
        "Rb",
        "Sr",
        "In",
        "Sn",
        "Sb",
        "Cs",
        "Ba",
        "Tl",
        "Pb",
        "Bi",
        "Po",
        "Fr",
        "Ra",
    )
)

# Terminal ligands that take the -1 when their bond to an over-valent centre is
# reduced, most electronegative first (they most naturally become the anion).
_TERMINAL_ANION_ELEMENTS = ("F", "O", "Cl", "N", "Br", "I", "S")

# One step down in bond order, relieving over-valence while keeping the atoms bonded.
# A single drops to UNSPECIFIED: a zero-valence bond (contributes 0, like ZERO) - the
# "zero-order" shed elsewhere in this module. UNSPECIFIED rather than ZERO because only it
# survives the biotite round-trip (ZERO has no biotite equivalent -> ANY + a per-bond
# warning; peppr.charge round-trips a sanitized mol via rdkit_interface.from_mol).
_REDUCED_BOND_ORDER = {
    Chem.BondType.TRIPLE: Chem.BondType.DOUBLE,
    Chem.BondType.DOUBLE: Chem.BondType.SINGLE,
    Chem.BondType.SINGLE: Chem.BondType.UNSPECIFIED,
}


def sanitize(mol: Chem.Mol, max_fix_iterations: int = 100) -> None:
    """
    Sanitize a molecule, repairing the problems that make RDKit's :func:`SanitizeMol`
    fail on otherwise-valid heavy-atom input.

    Parameters
    ----------
    mol : Chem.rdchem.Mol
        The molecule to sanitize, modified in place.
    max_fix_iterations : int, optional
        The maximum number of fix passes per stage.

    Notes
    -----
    Runs in stages, in the order RDKit resolves them - valence, then kekulization, then
    the rest - fixing and re-checking each before moving on:

    - ``AtomValenceException``: relieve over-valence via a zero-order bond (into a
      valence-unlimited metal or an over-valent partner), an "onium" charge (N/O), a
      terminal-ligand charge (hypervalent P/S), or a borate charge (B); metal centres and
      carborane cages are tolerated. See :func:`_fix_valence`.
    - ``KekulizeException``: kekulize a ring heteroatom (N/O/S) neutral-first - reduce its
      bond to a coordinated metal, or add a ring-N ``[nH]`` - and only as a last resort
      invent an onium ``+1`` (pyridinium/pyrylium/thiopyrylium). See :func:`_fix_kekule`.
    """
    # Avoid getting RDKit logs during checks for chemistry problems
    # If problems remain we rather want to receive them via an exception
    # from the final `SanitizeMol()` call
    with rdkit.rdBase.BlockLogs():
        # baseline clean-up first (standardises e.g. nitro groups and
        # organometallics); no corrections of ours are applied yet
        Chem.SanitizeMol(
            mol,
            sanitizeOps=SanitizeFlags.SANITIZE_CLEANUP
            | SanitizeFlags.SANITIZE_CLEANUP_ORGANOMETALLICS,
        )
        mol.UpdatePropertyCache(strict=False)

        # 1) valence problems must be resolved before anything else.
        _iteratively_fix(
            mol, SanitizeFlags.SANITIZE_PROPERTIES, _fix_valence, max_fix_iterations
        )
        try:
            Chem.SanitizeMol(mol, sanitizeOps=SanitizeFlags.SANITIZE_PROPERTIES)
        except Exception as exc:
            # tolerate over-valence RDKit still flags but no 2-centre-bond model can
            # fix - main-group metal centres and carborane cages, kept neutral (see
            # _is_tolerated_over_valence); anything else is a genuine over-valence
            # _fix_valence could not repair -> re-raise, noting which atoms
            offenders = [
                p.GetAtomIdx()
                for p in Chem.DetectChemistryProblems(
                    mol, sanitizeOps=SanitizeFlags.SANITIZE_PROPERTIES
                )
                if not _is_tolerated_over_valence(mol.GetAtomWithIdx(p.GetAtomIdx()))
            ]
            if offenders:
                exc.add_note(
                    f"peppr.sanitize could not fix over-valent atom(s) {offenders}"
                )
                raise

        # ring perception is required before kekulization can be attempted
        Chem.SanitizeMol(mol, sanitizeOps=SanitizeFlags.SANITIZE_SYMMRINGS)

        # 2) kekulization problems
        _iteratively_fix(
            mol, SanitizeFlags.SANITIZE_KEKULIZE, _fix_kekule, max_fix_iterations
        )
        # surface any remaining kekulization error here
        Chem.SanitizeMol(mol, sanitizeOps=SanitizeFlags.SANITIZE_KEKULIZE)

        # everything else (aromaticity, radicals, hybridisation, ...);
        # any issue that remains will be raised here
        Chem.SanitizeMol(
            mol,
            sanitizeOps=SanitizeFlags.SANITIZE_ALL ^ SanitizeFlags.SANITIZE_PROPERTIES,
        )


def _fix_valence(mol: Chem.Mol, problem: Exception) -> None:
    """
    Fix a single over-valence problem in place, applying the first repair that is
    legitimate for the atom.

    Parameters
    ----------
    mol : Chem.rdchem.Mol
        The molecule to fix.
    problem : Exception
        The problem to fix.

    Notes
    -----
    Strategies are attempted in order; the first that relieves the atom wins:

    0. Neutralise a spurious formal charge so the bond-based fixes see the real valence.
    1. A zero-order bond (:func:`_fix_valence_by_reduced_bond_order`) into a clean sink -
       a valence-unlimited metal (the FeMo-cofactor carbide bonded to 6 Fe) or an
       over-valent cage partner.
    2. Element-specific: N/O take an "onium" ``+1``; hypervalent P/S push a ``-1`` onto a
       terminal ligand; a non-cage tetravalent B becomes borate ``[B-]``; a carborane cage
       or main-group metal centre is tolerated neutral.

    An atom no fix relieves is a genuinely impossible valence (a divalent F, a pentavalent
    C), left flagged so :func:`sanitize` raises rather than papering over broken chemistry.
    """
    if problem.GetType() != "AtomValenceException":
        return
    at = mol.GetAtomWithIdx(problem.GetAtomIdx())

    # 0) drop a spurious formal charge first: biotite emits charged atoms illegal for
    #    their charge state (e.g. [C+2] carborane carbons, coordinated [Cl-]). Neutralise
    #    so the bond-based fixes below see the real valence; the element branches re-apply
    #    any legitimate charge (borate [B-], N/O onium, ...).
    if at.GetFormalCharge() != 0:
        at.SetFormalCharge(0)
        mol.UpdatePropertyCache(strict=False)

    # 1) decouple an over-valent bond as zero-order into a clean sink - a
    #    valence-unlimited metal (the FeMo-cofactor carbide bonded to 6 Fe) or a
    #    partner that is itself over-valent, e.g. a borane-cage neighbour (relieves
    #    both). No phantom H; direction-agnostic (no RWMol bond-flip).
    if _fix_valence_by_reduced_bond_order(mol, at):
        return

    # 2) element-specific charge, or tolerate an over-valence no 2-centre-bond model can
    #    express (a carborane cage or main-group metal centre - kept neutral, flagged,
    #    swallowed in `sanitize`). Each branch gates on the exact valence one over the
    #    element default, where an onium/borate charge is valid. Anything matching no branch
    #    is a genuinely impossible valence (divalent F, pentavalent C), left flagged so
    #    `sanitize` raises rather than papering it over with a phantom H.
    elem = at.GetSymbol()
    if elem == "O" and at.GetTotalValence() == 3:
        # oxonium: an over-coordinated O (one past its default valence 2) shares a lone
        # pair, relieved by a +1
        at.SetFormalCharge(1)
    elif elem == "N" and at.GetTotalValence() == 4:
        # ammonium/iminium: N one past its default valence 3, relieved by a +1
        at.SetFormalCharge(1)
        # a cumulated 1,3-dipole (azide/diazo) needs the compensating -1 on its terminal
        # N: CN=N=N -> CN=[N+]=[N-], not CN=[N+]=N
        _balance_dipole_terminal(at)
    elif elem in ("P", "S") and _fix_valence_by_ligand_charge(mol, at):
        # hypervalent P/S relieved by moving a -1 onto a terminal ligand: =O -> -O(-)
        # (phosphate, sulfonate/hypervalent sulfur), P-F -> P~[F-] (PF6). Nitro is
        # already handled by SANITIZE_CLEANUP. If no terminal ligand can take the charge,
        # this is False and the atom is left flagged to raise.
        pass
    elif elem == "B" and not _is_borane_cage_atom(at) and at.GetTotalValence() == 4:
        # a non-cage tetravalent boron is a borate anion [B-] (BF4-, B(OH)4-, R4B-) - the
        # N/O onium charge with the opposite sign.
        at.SetFormalCharge(-1)
    elif _is_tolerated_over_valence(at):
        # a carborane/borane cage vertex or a main-group metal centre (In, Be, ...): no
        # 2-centre-bond model expresses the hyper-coordination, so tolerate it neutral
        # (already neutralised in step 0, then flagged and swallowed in `sanitize` by the
        # same predicate).
        pass


def _fix_valence_by_reduced_bond_order(mol: Chem.Mol, atom: Chem.Atom) -> bool:
    """
    Relieve an over-valent atom by lowering a bond one order at a time into a sink that
    absorbs it (:data:`_REDUCED_BOND_ORDER`), one unit of excess per step. Return True if
    the atom ends up within valence.

    Notes
    -----
    A reduction debits *both* endpoints, so it is lossless only into a valence-unlimited
    metal (the FeMo carbide bonded to 6 Fe) or an over-valent cage partner (their shared
    bond's reduction relieves both). With no such sink the atom is left over-valent to
    raise later, deliberately avoiding a phantom H. Aromatic bonds are never touched.
    """
    if _has_no_valence_limit(atom.GetAtomicNum()):
        return False  # unenforced valence -> never flagged, nothing to relieve
    max_valence = max(Chem.GetPeriodicTable().GetValenceList(atom.GetAtomicNum()))
    excess = atom.GetTotalValence() - max_valence  # how far over its highest valence
    if excess <= 0:
        return False

    # spend the least-damaging bond first: a valence-unlimited metal (free sink), then an
    # over-valent cage partner, then an omitted main-group metal bond, then a covalent one.
    def _sink_rank(bond: Chem.Bond) -> int:
        partner = bond.GetOtherAtom(atom)
        neighbour = partner.GetAtomicNum()
        if _has_no_valence_limit(neighbour):
            return 0
        if (
            _is_over_valent(partner)
            and neighbour not in _ORGANOMETALLIC_OMITTED_ELEMENTS
        ):
            return 1
        if neighbour in _ORGANOMETALLIC_OMITTED_ELEMENTS:
            return 2
        return 3

    for bond in sorted(atom.GetBonds(), key=_sink_rank):
        if excess <= 0:
            break
        if bond.GetBondType() not in _REDUCED_BOND_ORDER:
            continue  # aromatic, already zero-order, dative, or otherwise not reducible
        # only reduce into a clean sink (absorbs it without a phantom H): a valence-
        # unlimited metal, or an over-valent partner in a cage cluster. The cage gate is
        # deliberate - an ordinary over-valent adjacency (pentavalent C next to over-valent
        # O) is broken input, left to raise, not silently zero-ordered like a cluster.
        partner = bond.GetOtherAtom(atom)
        clean_sink = _has_no_valence_limit(partner.GetAtomicNum()) or (
            _is_over_valent(partner)
            and (_is_borane_cage_atom(atom) or _is_borane_cage_atom(partner))
        )
        if not clean_sink:
            continue
        # drop one order (e.g. single -> zero), shedding exactly one valence unit per end
        bond.SetBondType(_REDUCED_BOND_ORDER[bond.GetBondType()])
        excess -= 1
        # refresh so a shared partner is not reduced past its own excess (no over-shedding)
        mol.UpdatePropertyCache(strict=False)

    mol.UpdatePropertyCache(strict=False)
    return bool(atom.GetTotalValence() <= max_valence)


def _fix_valence_by_ligand_charge(mol: Chem.Mol, atom: Chem.Atom) -> bool:
    """
    Relieve an over-valent centre by moving a ``-1`` onto its terminal ligands, keeping
    the centre itself neutral. Return True if resolved.

    Notes
    -----
    Per unit of excess, the bond to the most electronegative terminal (degree-1, neutral,
    H-free) ligand is lowered one order and that ligand takes the ``-1``: ``=O -> -O(-)``
    (phosphate) and ``P-F -> P~[F-]`` (PF6). Ligands are spent in electronegativity order
    (F > O > Cl > N > ...).
    """
    if _has_no_valence_limit(atom.GetAtomicNum()):
        return False
    max_valence = max(Chem.GetPeriodicTable().GetValenceList(atom.GetAtomicNum()))
    # amount by which the atom's valence exceeds its highest standard valence
    excess = atom.GetTotalValence() - max_valence
    if excess <= 0:
        return False

    rank = {sym: i for i, sym in enumerate(_TERMINAL_ANION_ELEMENTS)}
    # terminal, neutral, H-free ligands whose bond order can be reduced, most
    # electronegative first (they most naturally carry the -1)
    ligand_bonds = sorted(
        (
            bond
            for bond in atom.GetBonds()
            if bond.GetBondType() in _REDUCED_BOND_ORDER
            and bond.GetOtherAtom(atom).GetSymbol() in rank
            and bond.GetOtherAtom(atom).GetDegree() == 1
            and bond.GetOtherAtom(atom).GetFormalCharge() == 0
            and bond.GetOtherAtom(atom).GetTotalNumHs() == 0
        ),
        key=lambda bond: rank[bond.GetOtherAtom(atom).GetSymbol()],
    )
    for bond in ligand_bonds:
        if excess <= 0:
            break
        # lower the bond order by one and move the -1 onto the ligand
        terminal = bond.GetOtherAtom(atom)
        bond.SetBondType(_REDUCED_BOND_ORDER[bond.GetBondType()])
        excess -= 1
        terminal.SetFormalCharge(-1)

    mol.UpdatePropertyCache(strict=False)
    return bool(atom.GetTotalValence() <= max_valence)


def _balance_dipole_terminal(center: Chem.Atom) -> None:
    """
    Move the compensating ``-1`` onto the terminal N of a just-``+1``-charged cumulated
    1,3-dipole (azide/diazo): ``CN=N=N -> CN=[N+]=[N-]``.

    Notes
    -----
    Only a terminal (degree-1), neutral N double-bonded to ``center`` is charged, leaving
    ammonium/iminium centres untouched (RDKit would otherwise neutral-fill the divalent
    terminal N with an implicit H, giving ``CN=[N+]=N``).
    """
    for bond in center.GetBonds():
        if bond.GetBondType() != Chem.BondType.DOUBLE:
            continue
        terminal = bond.GetOtherAtom(center)
        if (
            terminal.GetSymbol() == "N"
            and terminal.GetDegree() == 1
            and terminal.GetFormalCharge() == 0
        ):
            terminal.SetFormalCharge(-1)
            return


def _has_no_valence_limit(atomic_num: int) -> bool:
    """
    True for elements whose valence RDKit does not enforce (transition and
    inner-transition metals). RDKit signals this with a ``[-1]`` sentinel in the
    valence list (not an empty list), so test for that sentinel explicitly.
    """
    valences = Chem.GetPeriodicTable().GetValenceList(atomic_num)
    return not valences or max(valences) < 0


def _is_over_valent(atom: Chem.Atom) -> bool:
    """
    True if the atom exceeds the larger of its neutral valence and its formally-charged
    (isoelectronic) valence - i.e. more bonds than any legitimate reading allows.

    Notes
    -----
    Taking the max reads a charge in whichever direction is legitimate:

    - a charge that *raises* the ceiling is honoured, so a ``[N+]`` (isoelectronic with C)
      or a ``[B-]`` (borate) may reach valence 4 without being over-valent (e.g. a BODIPY
      ``[N+]...[B-]``);
    - a spurious charge that would *lower* it is ignored, so biotite's ``[C+2]`` on a
      carborane carbon (isoelectronic with Be) is still judged by neutral carbon's 4 - it
      is dropped later anyway (step 0 of :func:`_fix_valence`).

    Valence-unlimited metals are never over-valent. Requires an up-to-date property cache.
    """
    if _has_no_valence_limit(atom.GetAtomicNum()):
        return False
    periodic_table = Chem.GetPeriodicTable()
    neutral_valence = max(periodic_table.GetValenceList(atom.GetAtomicNum()))
    # isoelectronic neutral element (same electron count): subtracting the charge shifts
    # the atomic number - left for a cation ([N+]->C), right for an anion ([B-]->C)
    isoelectronic = atom.GetAtomicNum() - atom.GetFormalCharge()
    charged_valence = (
        max(periodic_table.GetValenceList(isoelectronic)) if isoelectronic >= 1 else 0
    )
    return bool(atom.GetExplicitValence() > max(neutral_valence, charged_valence))


def _is_borane_cage_atom(atom: Chem.Atom) -> bool:
    """
    True if the atom is a borane / carborane cage vertex (bonded to >=3 borons).

    Notes
    -----
    Cage 3-centre-2-electron bonding is not expressible with 2-centre bonds: a cage carbon
    is dropped to a legal valence by the zero-order shed, while a cage boron has no clean
    fix and is tolerated neutral. A lone borate boron (bonded to F/O/C, no boron
    neighbours) is not a cage vertex and is charged ``[B-]`` instead. Boron-specific: the
    only element forming such cages in the PDB CCD.
    """
    return bool(
        sum(neighbour.GetSymbol() == "B" for neighbour in atom.GetNeighbors()) >= 3
    )


def _is_tolerated_over_valence(atom: Chem.Atom) -> bool:
    """
    True for over-valence peppr accepts as flagged-but-faithful: an omitted main-group
    metal centre (kept neutral) or a borane / carborane cage vertex
    (:func:`_is_borane_cage_atom`).

    Notes
    -----
    Both are hyper-coordinations no 2-centre-bond model can express. Everything else still
    over-valent after the fixes is a genuine, unrepairable error.
    """
    return (
        atom.GetAtomicNum() in _ORGANOMETALLIC_OMITTED_ELEMENTS
        or _is_borane_cage_atom(atom)
    )


def _fix_kekule(mol: Chem.Mol, problem: Exception) -> None:
    """
    Fix a single kekulization problem in place by pinning down the ambiguous protonation,
    charge, or metal-coordination state of a ring heteroatom (N/O/S).

    Parameters
    ----------
    mol : Chem.rdchem.Mol
        The molecule to fix.
    problem : Exception
        The problem to fix.

    Notes
    -----
    Ring perception (``SANITIZE_SYMMRINGS``) must have run first. The fixes are tried
    *neutral-first*, so a cation is only ever invented as a last resort:

    1. a heteroatom *coordinating a metal* (a pyridine/furan/thiophene donor) has an extra
       bond to the metal blocking kekulization; reduce that bond to zero-order, freeing the
       ring with the atom left neutral (and no directional dative bond, which ``sanitize``
       cannot orient editing in place);
    2. a two-neighbour ring N given an *explicit* H is read as a charge (``[nH]`` -> a
       charged ``[nH+]``), honouring an over-protonated input; heavy-atom input has only
       implicit H, so this never fires there;
    3. otherwise a two-neighbour ring N may just need an added implicit H - the pyrrole
       ``[nH]`` tautomer; nearly every unkekulizable heavy-atom ring (thiophenes, lactams,
       ...) is neutral and resolves here;
    4. only a saturated heteroatom that *no* neutral kekule can satisfy is a genuine onium,
       forced to carry a ring double bond, and takes a ``+1``: pyridinium (N, 3 neighbours),
       pyrylium (O, 2), or thiopyrylium (S, 2).
    """

    def _kekulized() -> bool:
        remaining = Chem.DetectChemistryProblems(
            mol, sanitizeOps=SanitizeFlags.SANITIZE_KEKULIZE
        )
        return not [p for p in remaining if _is_same_problem(p, problem)]

    if problem.GetType() != "KekulizeException":
        # unclear if it can be fixed
        return
    # gather the whole ring system around the flagged atoms
    problem_indices = set(problem.GetAtomIndices())
    ring_info = mol.GetRingInfo()
    for ring_indices in ring_info.AtomRings():
        if problem_indices.intersection(ring_indices):
            problem_indices |= set(ring_indices)

    # heavy-neighbour count at which a neutral aromatic heteroatom is saturated (lone pair
    # donated to the ring); past that it kekulizes only as an onium via a +1. N takes one
    # bond more than O/S.
    onium_saturation = {"N": 3, "O": 2, "S": 2}

    # candidate ring heteroatoms (N/O/S), grouped by the degree of freedom each carries
    two_neighbour_ns = []
    original_nH_state = []
    original_charge_state = []
    onium_atoms = []
    metal_coordinating_atoms = []
    for atidx in problem_indices:
        at = mol.GetAtomWithIdx(atidx)
        sym = at.GetSymbol()
        if sym not in onium_saturation:  # only ring N/O/S carry a kekulization DOF
            continue
        # aromatic heteroatom bonded to a metal: a neutral dative donor (coordinating
        # pyridine/furan/thiophene), relieved by reducing that bond, not by a charge
        if at.GetIsAromatic() and any(
            _is_metal(nbr.GetAtomicNum()) for nbr in at.GetNeighbors()
        ):
            metal_coordinating_atoms.append(at)
            continue
        num_heavy_neighbors = at.GetTotalDegree() - at.GetTotalNumHs()
        # a two-neighbour ring N is protonation-ambiguous (pyrrole [nH] vs pyridine [n])
        if sym == "N" and num_heavy_neighbors == 2:
            two_neighbour_ns.append(at)
            original_nH_state.append(at.GetNumExplicitHs())
            original_charge_state.append(at.GetFormalCharge())
        # a saturated neutral aromatic heteroatom is an onium (pyridinium/pyrylium/
        # thiopyrylium), needing a +1
        elif (
            at.GetIsAromatic()
            and at.GetTotalNumHs() == 0
            and at.GetFormalCharge() == 0
            and num_heavy_neighbors == onium_saturation[sym]
        ):
            onium_atoms.append(at)

    if not two_neighbour_ns and not onium_atoms and not metal_coordinating_atoms:
        # no candidate heteroatoms to fix the problem
        return

    # 1) metal-coordinating heteroatom: drop its bond to the metal to zero-order (the valence
    #    fix's relief into a metal sink, :data:`_REDUCED_BOND_ORDER`). That extra bond, not a
    #    missing charge/H, blocks kekulization; reducing it frees the ring, atom left neutral.
    reduced_metal_bonds = []  # (bond, original type) - restored if nothing resolves below
    for at in metal_coordinating_atoms:
        for nbr in at.GetNeighbors():
            if not _is_metal(nbr.GetAtomicNum()):
                continue
            bond = mol.GetBondBetweenAtoms(at.GetIdx(), nbr.GetIdx())
            if bond.GetBondType() in _REDUCED_BOND_ORDER:
                reduced_metal_bonds.append((bond, bond.GetBondType()))
                bond.SetBondType(_REDUCED_BOND_ORDER[bond.GetBondType()])
    if reduced_metal_bonds:
        mol.UpdatePropertyCache(strict=False)
        if _kekulized():
            return

    # 2) two-neighbour N given an *explicit* H: read that H as a charge ([nH] -> [nH+]),
    #    honouring an over-protonated input (an uncharged [nH]...[nH] imidazole is really the
    #    imidazolium) rather than stripping it. Only fires when the input supplied an explicit
    #    H, so the heavy-atom CCD path (implicit H) never reaches it - no CCD over-charging.
    for ni, nH, charge in zip(
        two_neighbour_ns, original_nH_state, original_charge_state
    ):
        if nH > 0 and charge == 0:
            ni.SetFormalCharge(nH)
            if _kekulized():
                return
            ni.SetFormalCharge(0)  # not resolved - try other candidates
    mol.UpdatePropertyCache(strict=False)

    # 3) two-neighbour N: search [nH] on/off - the *neutral* resolution (pyrrole [nH] donates
    #    its lone pair vs pyridine [n]). Tried before any charge: a ring that kekulizes with an
    #    added H needs no invented cation, so this must precede the onium (step 4). Nearly all
    #    unkekulizable heavy-atom rings (thiophenes, lactams, ...) resolve neutrally here.
    for combo in product([0, 1], repeat=len(two_neighbour_ns)):
        for ni, nH_explicit in zip(two_neighbour_ns, combo):
            ni.SetNumExplicitHs(nH_explicit)
        if _kekulized():
            return
    for ni, nH_explicit in zip(two_neighbour_ns, original_nH_state):
        ni.SetNumExplicitHs(nH_explicit)  # restore before the charge route
    mol.UpdatePropertyCache(strict=False)

    # 4) onium (last resort): only a genuine aromatic cation reaches here - a saturated ring
    #    heteroatom that no neutral kekule can satisfy, so it is forced to carry a ring double
    #    bond (pyridinium/pyrylium/thiazolium). Charge the fewest atoms, N before O/S so a
    #    mixed ring charges the right one (a thiazolium is an N-ylide, not an S+).
    onium_atoms.sort(key=lambda at: 0 if at.GetSymbol() == "N" else 1)
    for size in range(1, len(onium_atoms) + 1):
        for subset in combinations(onium_atoms, size):
            for at in subset:
                at.SetFormalCharge(1)
            if _kekulized():
                return
            for at in subset:
                at.SetFormalCharge(0)

    # nothing resolved it - roll back the metal reductions (the [nH] search already restored
    # its Hs above) so a failed attempt leaves the caller's in-place molecule untouched.
    for bond, original_type in reduced_metal_bonds:
        bond.SetBondType(original_type)
    mol.UpdatePropertyCache(strict=False)


def _is_metal(atomic_num: int) -> bool:
    """
    True for a metal peppr treats as a coordination centre: a transition / inner-transition
    metal (valence unenforced, :func:`_has_no_valence_limit`) or an omitted main-group metal
    (:data:`_ORGANOMETALLIC_OMITTED_ELEMENTS`). A ring heteroatom bonded to one is a neutral
    dative donor, not an onium.
    """
    return (
        _has_no_valence_limit(atomic_num)
        or atomic_num in _ORGANOMETALLIC_OMITTED_ELEMENTS
    )


def _iteratively_fix(
    mol: Chem.Mol,
    sanitize_op: SanitizeFlags,
    fix: Callable[[Chem.Mol, Exception], None],
    max_fix_iterations: int,
) -> None:
    """
    Repeatedly detect and fix the problems of a single sanitation stage, until every
    problem is either fixed or known to be unfixable here.

    Parameters
    ----------
    mol : Chem.rdchem.Mol
        The molecule to fix, modified in place.
    sanitize_op : SanitizeFlags
        The sanitation operation whose problems are detected and fixed.
    fix : Callable[[Chem.Mol, Exception], None]
        Function applying a fix for a single detected problem.
    max_fix_iterations : int
        The maximum number of iterations to fix problems.

    Notes
    -----
    RDKit surfaces problems iteratively - fixing one may resolve or shift others. Each pass
    re-detects, marks any problem still present after its fix was attempted as *stuck*
    (unfixable here, e.g. a tolerated metal centre), and acts only on the rest. Skipping
    stuck ones - rather than aborting the stage on the first recurrence - lets the rest of
    a cluster resolve completely.
    """
    stuck: list[Exception] = []
    attempted: list[Exception] = []
    for _ in range(max_fix_iterations):
        problems = Chem.DetectChemistryProblems(mol, sanitizeOps=sanitize_op)
        # a problem still present after its fix was attempted last pass cannot be
        # fixed here - record it so it is not retried and does not stall the stage
        for problem in problems:
            if any(_is_same_problem(problem, a) for a in attempted) and not any(
                _is_same_problem(problem, s) for s in stuck
            ):
                stuck.append(problem)
        # act only on problems not yet addressed (and not given up on)
        pending = [
            problem
            for problem in problems
            if not any(_is_same_problem(problem, s) for s in stuck)
        ]
        if not pending:
            # everything is either fixed or known to be unfixable here
            break
        for problem in pending:
            fix(mol, problem)
        attempted = pending


def _is_same_problem(problem1: Exception, problem2: Exception) -> bool:
    """
    Check if two exceptions related to chemistry problems are the same.

    Parameters
    ----------
    problem1, problem2 : Exception
        The problems to compare.

    Returns
    -------
    same : bool
        True if the two problems are the same.
    """
    problem_type = problem1.GetType()
    if problem_type != problem2.GetType():
        return False
    elif problem_type == "AtomValenceException":
        if problem1.GetAtomIdx() != problem2.GetAtomIdx():
            return False
    elif problem_type == "KekulizeException":
        if (
            len(
                set(problem1.GetAtomIndices()).intersection(
                    set(problem2.GetAtomIndices())
                )
            )
            == 0
        ):
            return False
    return True
