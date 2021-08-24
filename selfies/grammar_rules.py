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
from selfies.utils.smiles_utils import parse_bond_smiles


def parse_atom_selfies(symbol: str) -> Optional[Tuple[Any, Atom]]:
    try:
        bond_info, atom_fac = _PARSE_ATOM_CACHE[symbol]
    except KeyError:
        bond_info, atom_fac = _parse_atom_selfies_no_cache(symbol)
        _PARSE_ATOM_CACHE[symbol] = (bond_info, atom_fac)

    atom = atom_fac()
    if atom.bonding_capacity < 0:
        return None  # too many Hs (e.g. [CH9]
    return bond_info, atom


def parse_branch_selfies(symbol: str) -> Optional[Tuple[int, int]]:
    try:
        return _PARSE_BRANCH_CACHE[symbol]
    except KeyError:
        return None


def parse_ring_selfies(symbol: str) -> Optional[Tuple[Any, int]]:
    try:
        return _PARSE_RING_CACHE[symbol]
    except KeyError:
        return None


def next_atom_state(
        bond_order: int, bond_cap: int, state: int
) -> Tuple[int, int]:
    if state == -1:
        return 1, bond_cap
    bond_order = min(bond_order, state, bond_cap)
    next_state = bond_cap - bond_order
    return bond_order, next_state


def next_branch_state(branch_type: int, state: int) -> Tuple[int, int]:
    assert 1 <= branch_type <= 3

    if 2 <= state:
        branch_init_state = min(state - 1, branch_type)
        next_state = state - branch_init_state
        return branch_init_state, next_state
    else:
        return -1, state


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


def _parse_atom_selfies_no_cache(symbol):
    m = SELFIES_ATOM_PATTERN.match(symbol)
    if m is None:
        return None
    bond_char, isotope, element, chirality, h_count, charge = m.groups()

    if symbol[1 + len(bond_char):-1] in ORGANIC_SUBSET:
        atom_fac = functools.partial(Atom, element=element, is_aromatic=False)
        return parse_bond_smiles(bond_char), atom_fac

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

    return parse_bond_smiles(bond_char), atom_fac


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
        cache[symbol] = _parse_atom_selfies_no_cache(symbol)
    return cache


def _build_branch_cache():
    cache = dict()
    for L in range(1, 4):
        for bond_char in ["", "=", "#"]:
            symbol = "[{}Branch{}]".format(bond_char, L)
            cache[symbol] = (parse_bond_smiles(bond_char)[0], L)
    return cache


def _build_ring_cache():
    cache = dict()
    for L in range(1, 4):
        # [RingL], [=RingL], [#RingL]
        for bond_char in ["", "=", "#"]:
            symbol = "[{}Ring{}]".format(bond_char, L)
            bond_info = parse_bond_smiles(bond_char)
            cache[symbol] = ((bond_info, bond_info), L)

        # [-/RingL], [\/RingL], [\-RingL], ...
        for lchar, rchar in itertools.product(["-", "/", "\\"], repeat=2):
            if lchar == rchar == "-":
                continue
            symbol = "[{}{}Ring{}]".format(lchar, rchar, L)
            lbond_info = parse_bond_smiles(lchar)
            rbond_info = parse_bond_smiles(rchar)
            cache[symbol] = ((lbond_info, rbond_info), L)
    return cache


_PARSE_ATOM_CACHE = _build_atom_cache()

_PARSE_BRANCH_CACHE = _build_branch_cache()

_PARSE_RING_CACHE = _build_ring_cache()
