# mypy: disable-error-code="attr-defined"

"""Vendored copy of peppr's ``sanitize`` with local over-valence fixes.

.. warning::

    TEMPORARY. This is a local, hand-maintained copy of the modified
    ``sanitize`` from *peppr-internal*. It carries boron-cage / main-group
    over-valence tolerance and related valence fixes that are **not yet in a
    released** ``peppr`` and are not pip-installable from the internal repo.

    TODO(peppr): once these changes land in a published ``peppr`` release,
    DELETE this module and revert every import of the form
    ``from plinder.core.utils.sanitize import sanitize as peppr_sanitize``
    back to ``from peppr import sanitize as peppr_sanitize``
    (grep the tree for ``peppr_sanitize`` to find all call sites).
"""

__all__ = ["sanitize"]

from collections.abc import Callable
from itertools import product

import rdkit
import rdkit.Chem.AllChem as Chem
from rdkit.Chem.rdmolops import SanitizeFlags

# Main-group metals/metalloids that RDKit's SANITIZE_CLEANUP_ORGANOMETALLICS step
# leaves alone (it dativises only transition-metal coordination). An over-valent
# *centre* of one of these is tolerated neutral (see :func:`_fix_valence` and
# :func:`_is_tolerated_over_valence`), and a bond to one is a preferred zero-order
# sink over a genuine covalent bond (see :func:`_fix_valence_by_zero_order`).
# Atomic numbers: B, Li, Be, Na, Mg, Al, K, Ca, Ga, Ge, As, Rb, Sr, In, Sn, Sb, Cs,
# Ba, Tl, Pb, Bi, Po, Fr, Ra.
_ORGANOMETALLIC_OMITTED_ELEMENTS = frozenset(
    {
        5,
        3,
        4,
        11,
        12,
        13,
        19,
        20,
        31,
        32,
        33,
        37,
        38,
        49,
        50,
        51,
        55,
        56,
        81,
        82,
        83,
        84,
        87,
        88,
    }
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
    True if the atom has more bonds than *any* reasonable reading of it allows - the
    larger of its neutral valence and its formally-charged (isoelectronic) valence.
    Taking the max reads a charge in whichever direction is legitimate:

    - a charge that *raises* the ceiling is honoured, so a ``[N+]`` (isoelectronic with C)
      or a ``[B-]`` (borate) may reach valence 4 without being over-valent - and so is not
      over-shed into a phantom H (e.g. a BODIPY ``[N+]...[B-]``);
    - a spurious input charge that would *lower* the ceiling is ignored, so biotite's
      ``[C+2]`` on a carborane carbon (isoelectronic with Be, max valence 2) is still
      judged by neutral carbon's 4 - a cage boron then stops shedding into it at valence 4
      instead of draining it to a phantom H. The ``[C+2]`` is dropped later anyway (step 0
      of :func:`_fix_valence`), so honouring the neutral ceiling only anticipates that.

    Valence-unlimited metals are never over-valent (RDKit does not enforce them).
    Requires an up-to-date property cache, as it reads the explicit valence.
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
    True if the atom is a borane / carborane cage vertex (bonded to >=3 borons). Cage
    3-centre-2-electron bonding is not expressible with 2-centre bonds: a cage carbon is
    dropped to a legal valence by the zero-order shed, while a cage boron has no clean
    fix and is tolerated neutral. A lone borate boron (bonded to F/O/C, no boron
    neighbours) is not a cage vertex and is charged to ``[B-]`` instead. Boron-specific:
    the only element forming such cages in the PDB CCD.
    """
    return bool(
        sum(neighbour.GetSymbol() == "B" for neighbour in atom.GetNeighbors()) >= 3
    )


def _is_tolerated_over_valence(atom: Chem.Atom) -> bool:
    """
    True for over-valence peppr accepts as flagged-but-faithful rather than a
    failure: an omitted main-group metal *centre* (kept neutral) or a borane /
    carborane cage vertex (:func:`_is_borane_cage_atom`) - both are
    hyper-coordinations no 2-centre-bond model can express. Everything else still
    over-valent after the fixes is a genuine, unrepairable error.
    """
    return (
        atom.GetAtomicNum() in _ORGANOMETALLIC_OMITTED_ELEMENTS
        or _is_borane_cage_atom(atom)
    )


def sanitize(mol: Chem.Mol, max_fix_iterations: int = 100) -> None:
    """
    Fix small issues with RDKit SanitizeMol and sanitize molecule.

    This is an alternative to using :func:`SanitizeMol()` directly, in cases it fails
    due to issues with the molecule which are fixed by this function.

    Sanitation is run in stages so problems are fixed and surfaced in the order
    RDKit itself resolves them: valence first, then kekulization, then the rest.

    Parameters
    ----------
    mol : Chem.rdchem.Mol
        The molecule to sanitize.
    max_fix_iterations : int, optional
        The maximum number of iterations to fix problems.

    Notes
    -----
    Deals with cases:

    - ``AtomValenceException``: relieve over-valence via a zero-order bond (into a
      valence-unlimited metal or an over-valent partner), a formal "onium" charge
      (N/O), a terminal-ligand charge (hypervalent P/S), or a borate charge (B);
      metal centres and carborane cages are tolerated. See :func:`_fix_valence`.
    - ``KekulizeException``: for ring N atoms with unspecified protonation
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
    Apply fixes for common valence problems flagged by RDKit SanitizeMol.

    Strategies are attempted in order; the first that relieves the atom wins:

    0. Neutralise a spurious formal charge so the bond-based fixes see the real
       valence.
    1. Zero-order bond (:func:`_fix_valence_by_zero_order`) into a clean sink - a
       valence-unlimited (transition) metal (the FeMo-cofactor carbide bonded to
       6 Fe) or an over-valent partner (a borane-cage neighbour; this is how a cage
       carbon is dropped to a legal valence).
    2. Element-specific: N/O take a positive "onium" charge; hypervalent P/S push a
       ``-1`` onto a terminal ligand; a non-cage tetravalent boron becomes borate
       ``[B-]``; a carborane cage vertex or a main-group metal centre (In, Be, ...)
       is tolerated neutral (flagged, swallowed in :func:`sanitize`).

    An atom that none of these fixes relieves is over-valent past its element's maximum
    with no legitimate representation - a genuinely impossible valence (a divalent F, a
    pentavalent C) - and is left flagged so :func:`sanitize` raises rather than papering
    over broken chemistry.

    Parameters
    ----------
    mol : Chem.rdchem.Mol
        The molecule to fix.
    problem : Exception
        The problem to fix.
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
    if _fix_valence_by_zero_order(mol, at):
        return

    # 2) element-specific charge, or tolerate an over-valence no 2-centre-bond model can
    #    express (a carborane cage vertex or a main-group metal centre - kept neutral and
    #    flagged, swallowed in `sanitize`)
    # Each branch applies the one legitimate fix for its element. An onium/borate charge is
    # only valid one unit over the element default (the atom is neutral after step 0, so
    # "one over the default" is exactly `_is_over_valent` by +1), so those branches gate on
    # that exact valence. Anything that matches no branch is over-valent past its element's
    # maximum with no legitimate fix - a genuinely impossible valence (a divalent F, a
    # pentavalent C) - and is deliberately left flagged so `sanitize` raises it rather than
    # papering over broken chemistry with a phantom H.
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


def _fix_valence_by_zero_order(mol: Chem.Mol, atom: Chem.Atom) -> bool:
    """
    Relieve an over-valent atom by lowering its bonds one order at a time (a single
    becomes a zero-order bond via :data:`_REDUCED_BOND_ORDER`), one unit of excess per
    step. Reducing a bond needs no donor/acceptor direction and can relieve the bond's
    *end* atom (no ``RWMol`` bond-flip), but it debits *both* endpoints - so it is
    lossless only into a sink that can absorb it: a valence-unlimited metal (e.g. the
    FeMo-cofactor carbide bonded to 6 Fe) or a partner that is itself over-valent
    (the shared bond's reduction relieves both).

    Only such clean sinks are used, so no phantom H is ever injected; an atom with no
    clean sink is left over-valent (and, if no other branch of :func:`_fix_valence`
    fixes it, surfaces as a raise rather than being papered over). The least-damaging
    bond is spent first (see :func:`_sink_rank`). Returns True if ``atom`` is now within
    its standard valence.
    """
    if _has_no_valence_limit(atom.GetAtomicNum()):
        # the atom itself is unenforced -> it is never flagged, nothing to relieve
        return False
    max_valence = max(Chem.GetPeriodicTable().GetValenceList(atom.GetAtomicNum()))
    # amount by which the atom's valence exceeds its highest standard valence
    excess = atom.GetTotalValence() - max_valence
    if excess <= 0:
        return False

    # spend the least-damaging bonds first, and relieve a hard-to-relieve partner
    # before an easy one: a valence-unlimited metal is a free sink; then an over-valent
    # *non-omitted* partner (a cage carbon - it strands if the borons satisfy each other
    # first); then an omitted main-group metal (coordinate) bond; a covalent bond last.
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
            continue  # already zero-order, dative, or otherwise not reducible
        # only reduce into a clean sink - one that absorbs the reduction without a phantom
        # H: a valence-unlimited metal (no default valence to fill) or an over-valent
        # partner in a borane cage (the shared bond's reduction relieves it too). The cage
        # gate is deliberate: an ordinary over-valent adjacency (e.g. a pentavalent C next
        # to an over-valent O) is broken input and is left to raise, not silently
        # zero-ordered the way a cluster is. Anything else would under-fill the partner
        # into a phantom H, so it is left alone.
        partner = bond.GetOtherAtom(atom)
        clean_sink = _has_no_valence_limit(partner.GetAtomicNum()) or (
            _is_over_valent(partner)
            and (_is_borane_cage_atom(atom) or _is_borane_cage_atom(partner))
        )
        if not clean_sink:
            continue
        # step the bond down one order (single -> zero-order), removing exactly one
        # unit of valence from both ends so a double/triple is not over-decoupled
        bond.SetBondType(_REDUCED_BOND_ORDER[bond.GetBondType()])
        excess -= 1
        # refresh so the next iteration sees the partner's reduced valence and will
        # not reduce a shared partner past its own excess (avoids over-shedding)
        mol.UpdatePropertyCache(strict=False)

    mol.UpdatePropertyCache(strict=False)
    return bool(atom.GetTotalValence() <= max_valence)


def _fix_valence_by_ligand_charge(mol: Chem.Mol, atom: Chem.Atom) -> bool:
    """
    Relieve an over-valent centre by moving charge onto its terminal ligands, keeping
    the centre neutral. Per unit of excess, the bond to the most electronegative
    terminal (degree-1, neutral, H-free) ligand is lowered one order and that ligand
    takes a ``-1``: ``=O -> -O(-)`` (phosphate) and ``P-F -> P~[F-]`` (PF6). Ligands are
    spent in electronegativity order (F > O > Cl > N > ...). Returns True if resolved.
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
    1,3-dipole (azide/diazo): ``CN=N=N -> CN=[N+]=[N-]``, not ``CN=[N+]=N`` (where RDKit
    would neutral-fill the divalent terminal N with an implicit H). Only a terminal
    (degree-1), neutral N double-bonded to ``center`` is charged, leaving
    ammonium/iminium centres untouched.
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


def _fix_kekule(mol: Chem.Mol, problem: Exception) -> None:
    """
    Apply fixes for common kekulization problems flagged by RDKit SanitizeMol.

    Targets ring nitrogens with unspecified protonation/charge (the usual cause
    of a ``KekulizeException``): first try interpreting an existing explicit H as
    a formal charge, then search combinations of ``[nH]`` assignments. Ring
    perception (``SANITIZE_SYMMRINGS``) must have run before this is called.

    Parameters
    ----------
    mol : Chem.rdchem.Mol
        The molecule to fix.
    problem : Exception
        The problem to fix.
    """
    if problem.GetType() != "KekulizeException":
        # unclear if it can be fixed
        return
    # get all atoms in the same ring as the problem atoms
    problem_indices = set(problem.GetAtomIndices())
    ring_info = mol.GetRingInfo()
    for ring_indices in ring_info.AtomRings():
        if problem_indices.intersection(ring_indices):
            problem_indices |= set(ring_indices)
    # if KekulizeException is related to ring N atoms - it could be due to
    # either incorrectly specified explicit protonation or charge
    original_nH_state = []
    original_charge_state = []
    possible_nHs = []
    for atidx in problem_indices:
        at = mol.GetAtomWithIdx(atidx)
        num_heavy_neighbors = at.GetTotalDegree() - at.GetTotalNumHs()
        if at.GetSymbol() == "N" and num_heavy_neighbors == 2:
            original_nH_state.append(at.GetNumExplicitHs())
            original_charge_state.append(at.GetFormalCharge())
            possible_nHs.append(at)

    if not len(possible_nHs):
        # no candidate nitrogens to fix the problem
        # note this early return is unnecessary as the loops below will
        # not be entered, but it makes the intention clearer
        return

    # first test if the issue is the formal charge due to explicit protonation
    for ni, nH, charge in zip(possible_nHs, original_nH_state, original_charge_state):
        if nH > 0 and charge == 0:
            # if we found a protonated nitrogen, we can fix the charge
            ni.SetFormalCharge(nH)
            # and check if the problem is resolved
            check_problems = Chem.DetectChemistryProblems(
                mol, sanitizeOps=SanitizeFlags.SANITIZE_KEKULIZE
            )
            if not [p for p in check_problems if _is_same_problem(p, problem)]:
                # problem is resolved - no need to try other combinations
                return
            else:
                # if the problem is not resolved - reset the charge to try other candidates
                ni.SetFormalCharge(charge)

    # alternatively, attempt to fix by setting one or more ring nitrogens as "[nH1]"
    # if there is one or more candidates - iterate through all the possible
    # combinations of explicit Hs for the candidate nitrogens
    for combo in product([0, 1], repeat=len(possible_nHs)):
        for ni, nH_explicit in zip(possible_nHs, combo):
            ni.SetNumExplicitHs(nH_explicit)
        # and check if the problem is resolved
        check_problems = Chem.DetectChemistryProblems(
            mol, sanitizeOps=SanitizeFlags.SANITIZE_KEKULIZE
        )
        if not [p for p in check_problems if _is_same_problem(p, problem)]:
            # problem is resolved - no need to try other combinations
            return
    # if all combos failed - reset explicit Hs back to the original values
    for ni, nH_explicit in zip(possible_nHs, original_nH_state):
        ni.SetNumExplicitHs(nH_explicit)


def _iteratively_fix(
    mol: Chem.Mol,
    sanitize_op: SanitizeFlags,
    fix: Callable[[Chem.Mol, Exception], None],
    max_fix_iterations: int,
) -> None:
    """
    Repeatedly detect and fix the problems of a single sanitation stage.

    RDKit surfaces problems iteratively - fixing one may resolve or shift others.
    Each pass re-detects, marks any problem that is still present after its fix
    was attempted as *stuck* (unfixable here, e.g. a tolerated metal centre), and
    fixes only the remaining, not-yet-addressed problems. Skipping the stuck ones
    - rather than aborting the whole stage on the first recurrence - lets the rest
    of a cluster be resolved completely (including problems newly shifted onto
    neighbours). Stops when every problem is either fixed or stuck, or
    ``max_fix_iterations`` is reached.

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
