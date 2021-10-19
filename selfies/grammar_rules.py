import functools
import itertools
import re
from typing import Any, List, Optional, Tuple

from selfies.constants import (
    ELEMENTS,
    INDEX_ALPHABET,
    INDEX_CODE,
    ORGANIC_SUBSET
)
from selfies.mol_graph import Atom
from selfies.utils.smiles_utils import smiles_to_bond


def process_atom_symbol(symbol: str) -> Optional[Tuple[Any, Atom]]:
    try:
        output = _PROCESS_ATOM_CACHE[symbol]
    except KeyError:
        output = _process_atom_selfies_no_cache(symbol)
        if output is None:
            return None
        _PROCESS_ATOM_CACHE[symbol] = output

    bond_info, atom_fac = output
    atom = atom_fac()
    if atom.bonding_capacity < 0:
        return None  # too many Hs (e.g. [CH9]
    return bond_info, atom


def process_branch_symbol(symbol: str) -> Optional[Tuple[int, int]]:
    try:
        return _PROCESS_BRANCH_CACHE[symbol]
    except KeyError:
        return None


def process_ring_symbol(symbol: str) -> Optional[Tuple[int, int, Any]]:
    try:
        return _PROCESS_RING_CACHE[symbol]
    except KeyError:
        return None


def next_atom_state(
        bond_order: int, bond_cap: int, state: int
) -> Tuple[int, Optional[int]]:
    if state == 0:
        bond_order = 0

    bond_order = min(bond_order, state, bond_cap)
    bonds_left = bond_cap - bond_order
    next_state = None if (bonds_left == 0) else bonds_left
    return bond_order, next_state


def next_branch_state(
        branch_type: int, state: int
) -> Tuple[int, Optional[int]]:
    assert 1 <= branch_type <= 3
    assert state > 1

    branch_init_state = min(state - 1, branch_type)
    next_state = state - branch_init_state
    return branch_init_state, next_state


def next_ring_state(
        ring_type: int, state: int
) -> Tuple[int, Optional[int]]:
    assert state > 0

    bond_order = min(ring_type, state)
    bonds_left = state - bond_order
    next_state = None if (bonds_left == 0) else bonds_left
    return bond_order, next_state


def get_index_from_selfies(*symbols: List[str]) -> int:
    index = 0
    for i, c in enumerate(reversed(symbols)):
        index += INDEX_CODE.get(c, 0) * (len(INDEX_CODE) ** i)
    return index


def get_selfies_from_index(index: int) -> List[str]:
    if index < 0:
        raise IndexError()
    elif index == 0:
        return [INDEX_ALPHABET[0]]

    symbols = []
    base = len(INDEX_ALPHABET)
    while index:
        symbols.append(INDEX_ALPHABET[index % base])
        index //= base
    return symbols[::-1]


# =============================================================================
# Caches (for computational speed)
# =============================================================================


SELFIES_ATOM_PATTERN = re.compile(
    r"^[\[]"  # opening square bracket [
    r"([=#/\\]?)"  # bond char
    r"(\d*)"  # isotope number (optional, e.g. 123, 26)
    r"([A-Z][a-z]?)"  # element symbol
    r"([@]{0,2})"  # chiral_tag (optional, only @ and @@ supported)
    r"((?:[H]\d)?)"  # H count (optional, e.g. H1, H3)
    r"((?:[+-][1-9]+)?)"  # charge (optional, e.g. +1)
    r"[]]$"  # closing square bracket ]
)


def _process_atom_selfies_no_cache(symbol):
    m = SELFIES_ATOM_PATTERN.match(symbol)
    if m is None:
        return None
    bond_char, isotope, element, chirality, h_count, charge = m.groups()

    if symbol[1 + len(bond_char):-1] in ORGANIC_SUBSET:
        atom_fac = functools.partial(Atom, element=element, is_aromatic=False)
        return smiles_to_bond(bond_char), atom_fac

    isotope = None if (isotope == "") else int(isotope)
    if element not in ELEMENTS:
        return None
    chirality = None if (chirality == "") else chirality

    s = h_count
    if s == "":
        h_count = 0
    else:
        h_count = int(s[1:])

    s = charge
    if s == "":
        charge = 0
    else:
        charge = int(s[1:])
        charge *= 1 if (s[0] == "+") else -1

    atom_fac = functools.partial(
        Atom,
        element=element,
        is_aromatic=False,
        isotope=isotope,
        chirality=chirality,
        h_count=h_count,
        charge=charge
    )

    return smiles_to_bond(bond_char), atom_fac


def _build_atom_cache():
    cache = dict()
    common_symbols = [
        "[#C+1]", "[#C-1]", "[#C]", "[#N+1]", "[#N]", "[#O+1]", "[#P+1]",
        "[#P-1]", "[#P]", "[#S+1]", "[#S-1]", "[#S]", "[=C+1]", "[=C-1]",
        "[=C]", "[=N+1]", "[=N-1]", "[=N]", "[=O+1]", "[=O]", "[=P+1]",
        "[=P-1]", "[=P]", "[=S+1]", "[=S-1]", "[=S]", "[Br]", "[C+1]", "[C-1]",
        "[C]", "[Cl]", "[F]", "[H]", "[I]", "[N+1]", "[N-1]", "[N]", "[O+1]",
        "[O-1]", "[O]", "[P+1]", "[P-1]", "[P]", "[S+1]", "[S-1]", "[S]"
    ]

    for symbol in common_symbols:
        cache[symbol] = _process_atom_selfies_no_cache(symbol)
    return cache


def _build_branch_cache():
    cache = dict()
    for L in range(1, 4):
        for bond_char in ["", "=", "#"]:
            symbol = "[{}Branch{}]".format(bond_char, L)
            cache[symbol] = (smiles_to_bond(bond_char)[0], L)
    return cache


def _build_ring_cache():
    cache = dict()
    for L in range(1, 4):
        # [RingL], [=RingL], [#RingL]
        for bond_char in ["", "=", "#"]:
            symbol = "[{}Ring{}]".format(bond_char, L)
            order, stereo = smiles_to_bond(bond_char)
            cache[symbol] = (order, L, (stereo, stereo))

        # [-/RingL], [\/RingL], [\-RingL], ...
        for lchar, rchar in itertools.product(["-", "/", "\\"], repeat=2):
            if lchar == rchar == "-":
                continue
            symbol = "[{}{}Ring{}]".format(lchar, rchar, L)
            order, lstereo = smiles_to_bond(lchar)
            order, rstereo = smiles_to_bond(rchar)
            cache[symbol] = (order, L, (lstereo, rstereo))
    return cache


_PROCESS_ATOM_CACHE = _build_atom_cache()

_PROCESS_BRANCH_CACHE = _build_branch_cache()

_PROCESS_RING_CACHE = _build_ring_cache()
