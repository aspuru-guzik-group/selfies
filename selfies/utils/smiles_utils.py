import enum
import re
from collections import deque
from typing import Iterator, Optional, Tuple, Union

from selfies.constants import AROMATIC_SUBSET, ELEMENTS, ORGANIC_SUBSET
from selfies.exceptions import SMILESParserError
from selfies.mol_graph import Atom, DirectedBond, MolecularDFSTree

SMILES_BRACKETED_ATOM_PATTERN = re.compile(
    r"^[\[]"  # opening square bracket [
    r"(\d*)"  # isotope number (optional, e.g. 123, 26)
    r"([A-Za-z][a-z]?)"  # element symbol
    r"([@]{0,2})"  # chiral_tag (optional, only @ and @@ supported)
    r"((?:[H]\d?)?)"  # H count (optional, e.g. H, H0, H3)
    r"((?:[+]+|[-]+|[+-]\d+)?)"  # charge (optional, e.g. ---, +1, ++)
    r"((?:[:]\d+)?)"  # atom class (optional, e.g. :12, :1)
    r"[]]$"  # closing square bracket ]
)

SMILES_BOND_ORDERS = {"-": 1, "/": 1, "\\": 1, ":": 1.5, "=": 2, "#": 3}
SMILES_STEREO_BONDS = {"/", "\\"}


class SMILESTokenTypes(enum.Enum):
    ATOM = 0
    BRANCH = 1
    RING = 2
    DOT = 3


class SMILESToken:

    def __init__(
            self,
            bond_idx: Optional[int],
            start_idx: int, end_idx: int, token_type: SMILESTokenTypes
    ):
        self.bond_idx = bond_idx
        self.start_idx = start_idx
        self.end_idx = end_idx
        self.token_type = token_type

    def extract_bond_char(self, smiles):
        return None if (self.bond_idx is None) else smiles[self.bond_idx]

    def extract_symbol(self, smiles):
        return smiles[self.start_idx:self.end_idx]


def tokenize_smiles(smiles: str) -> Iterator[SMILESToken]:
    i = 0
    while i < len(smiles):

        if smiles[i] == ".":
            yield SMILESToken(None, i, i + 1, SMILESTokenTypes.DOT)
            i += 1
            continue

        if smiles[i] in SMILES_BOND_ORDERS:
            bond_idx = i
            i += 1
        else:
            bond_idx = None

        if i == len(smiles):
            raise SMILESParserError(smiles, "hanging bond", i)

        elif smiles[i].isalpha():  # organic subset elements
            if smiles[i: i + 2] in ("Br", "Cl"):  # two letter elements
                token = SMILESToken(bond_idx, i, i + 2, SMILESTokenTypes.ATOM)
            else:  # one letter elements (e.g. C, N, ...)
                token = SMILESToken(bond_idx, i, i + 1, SMILESTokenTypes.ATOM)

        elif smiles[i] == "[":  # atoms encased in brackets (e.g. [NH])
            r_idx = smiles.find("]", i + 1)
            if r_idx == -1:
                raise SMILESParserError(smiles, "hanging bracket [", i)
            token = SMILESToken(bond_idx, i, r_idx + 1, SMILESTokenTypes.ATOM)

        elif smiles[i] in ("(", ")"):  # open and closed branch brackets
            token = SMILESToken(None, i, i + 1, SMILESTokenTypes.BRANCH)

        elif smiles[i].isdigit():  # one-digit ring number
            token = SMILESToken(bond_idx, i, i + 1, SMILESTokenTypes.RING)

        elif smiles[i] == "%":  # two-digit ring number (e.g. %12)
            rnum = smiles[i + 1: i + 3]
            if not (rnum.isnumeric() and len(rnum) == 2):
                err_msg = "invalid ring number '%{}'".format(rnum)
                raise SMILESParserError(smiles, err_msg, i)
            token = SMILESToken(bond_idx, i, i + 3, SMILESTokenTypes.RING)

        else:
            err_msg = "unrecognized symbol '{}'".format(smiles[i])
            raise SMILESParserError(smiles, err_msg, i)

        yield token
        i = token.end_idx


# =============================================================================
# SMILES -> Atom, Graph, etc.
# =============================================================================


def smiles_to_atom(atom_symbol: str) -> Optional[Atom]:
    if atom_symbol[0] == "[" and atom_symbol[-1] == "]":
        pass  # continue below
    elif atom_symbol in ORGANIC_SUBSET:
        return Atom(atom_symbol, False)
    elif atom_symbol in AROMATIC_SUBSET:
        return Atom(atom_symbol.capitalize(), True)
    else:
        return None

    m = SMILES_BRACKETED_ATOM_PATTERN.match(atom_symbol)
    if m is None:
        return None
    isotope, element, chirality, h_count, charge, _ = m.groups()

    isotope = None if (isotope == "") else int(isotope)
    is_aromatic = element.islower() and (element in AROMATIC_SUBSET)
    element = element.capitalize()
    if element not in ELEMENTS:
        return None
    chirality = None if (chirality == "") else chirality

    s = h_count
    if s == "":
        h_count = 0
    else:
        s = s[1:]  # HXXX -> XXX
        h_count = 1 if (s == "") else int(s)

    s = charge
    if s == "":
        charge = 0
    else:
        if s[-1].isdigit():  # (+/-)XXX
            charge = int(s[1:])
        else:  # +++... or ---....
            charge = len(s)
        charge *= 1 if s[0] == "+" else -1

    return Atom(
        element=element,
        is_aromatic=is_aromatic,
        isotope=isotope,
        chirality=chirality,
        h_count=h_count,
        charge=charge
    )


def smiles_to_bond(
        bond_char: Optional[str]
) -> Tuple[Union[int, float], Optional[str]]:
    order = SMILES_BOND_ORDERS.get(bond_char, 1)
    stereo = bond_char if (bond_char in SMILES_STEREO_BONDS) else None
    return order, stereo


def smiles_to_mol(smiles: str) -> MolecularDFSTree:
    if smiles == "":
        raise SMILESParserError(smiles, "empty SMILES", 0)

    mol = MolecularDFSTree()
    tokens = deque(tokenize_smiles(smiles))
    while tokens:
        _derive_mol_from_tokens(mol, smiles, tokens)
    return mol


def _derive_mol_from_tokens(mol, smiles, tokens):
    tok = None
    prev_stack = deque()
    branch_stack = deque()
    ring_log = dict()
    chain_start = True

    prev_stack.append(tok)
    while tokens:
        tok = tokens.popleft()
        bond_char = tok.extract_bond_char(smiles)
        symbol, symbol_type = tok.extract_symbol(smiles), tok.token_type
        prev_atom = prev_stack[-1]

        if symbol_type == SMILESTokenTypes.DOT:
            break

        elif symbol_type == SMILESTokenTypes.ATOM:
            curr = smiles_to_atom(symbol)
            if curr is None:
                err_msg = "invalid atom symbol '{}'".format(symbol)
                raise SMILESParserError(smiles, err_msg, tok.start_idx)

            curr = _attach_atom(mol, bond_char, curr, prev_atom)
            prev_stack.pop()
            prev_stack.append(curr)
            chain_start = False

        elif chain_start:
            err_msg = "SMILES chain begins with non-atom"
            raise SMILESParserError(smiles, err_msg, tok.start_idx)

        elif symbol_type == SMILESTokenTypes.BRANCH:
            if symbol == "(":
                branch_stack.append(tok)
                prev_stack.append(prev_atom)
                chain_start = True
            else:
                if not branch_stack:
                    err_msg = "hanging ')' bracket"
                    raise SMILESParserError(smiles, err_msg, tok.start_idx)
                branch_stack.pop()
                prev_stack.pop()

        elif symbol_type == SMILESTokenTypes.RING:
            if symbol not in ring_log:
                lpos = mol.add_placeholder_bond(src=prev_atom.index)
                ring_log[symbol] = (tok, prev_atom, lpos)
            else:
                ltoken, latom, lpos = ring_log.pop(symbol)
                _make_ring_bonds(
                    mol=mol, smiles=smiles,
                    ltoken=ltoken, latom=latom, lpos=lpos,
                    rtoken=tok, ratom=prev_atom
                )

        else:
            # should not happen
            raise Exception("invalid symbol type")

    if len(mol) == 0:
        err_idx = (len(smiles) if (tok is None) else tok.start_idx) - 1
        raise SMILESParserError(smiles, "empty SMILES fragment", err_idx)

    if branch_stack:
        err_idx = branch_stack[-1].start_idx
        raise SMILESParserError(smiles, "hanging '(' bracket", err_idx)

    if ring_log:
        rnum, (tok, _, _) = list(ring_log.items())[-1]
        err_msg = "hanging ring number '{}'".format(rnum)
        raise SMILESParserError(smiles, err_msg, tok.start_idx)


def _attach_atom(mol, bond_char, atom, prev_atom):
    is_root = (prev_atom is None)
    mol.add_atom(atom, mark_root=is_root)

    if not is_root:
        src, dst = prev_atom.index, atom.index
        order, stereo = smiles_to_bond(bond_char)
        if prev_atom.is_aromatic and atom.is_aromatic and (bond_char is None):
            order = 1.5
        mol.add_bond(src=src, dst=dst, order=order, stereo=stereo)
    return atom


def _make_ring_bonds(mol, smiles, ltoken, latom, lpos, rtoken, ratom):
    if mol.has_bond(latom.index, ratom.index):
        err_msg = "ring bond specified between already-bonded atoms"
        raise SMILESParserError(smiles, err_msg, ltoken.start_idx)

    lbond_char = ltoken.extract_bond_char(smiles)
    rbond_char = rtoken.extract_bond_char(smiles)

    # checking that ring bonds match
    bonds = (lbond_char, rbond_char)
    if bonds[0] is None:  # swap tuple elements
        bonds = (bonds[1], bonds[0])

    if ((bonds[0] == bonds[1])
            or (bonds[1] is None)
            or all(x in SMILES_STEREO_BONDS for x in bonds)):
        pass
    else:
        err_msg = "mismatched ring bonds"
        raise SMILESParserError(smiles, err_msg, ltoken.start_idx)

    lorder, lstereo = smiles_to_bond(lbond_char)
    rorder, rstereo = smiles_to_bond(rbond_char)
    if latom.is_aromatic and ratom.is_aromatic and (bonds == (None, None)):
        lorder = rorder = 1.5

    mol.add_ring_bond(
        a=latom.index, a_stereo=lstereo, a_pos=lpos,
        b=ratom.index, b_stereo=rstereo,
        order=max(lorder, rorder)
    )


# =============================================================================
# SMILES <- Atom, Graph, etc.
# =============================================================================


def atom_to_smiles(atom: Atom, brackets: bool = True) -> str:
    assert not atom.is_aromatic

    specs = (atom.isotope, atom.chirality, atom.h_count, atom.charge)
    if specs == (None, None, None, 0):
        return atom.element
    else:
        builder = []
        if brackets:
            builder.append("[")
        if atom.isotope is not None:
            builder.append(str(atom.isotope))
        builder.append(atom.element)
        if atom.chirality is not None:
            builder.append(atom.chirality)
        if atom.h_count != 0:
            builder.append("H")
            builder.append(str(atom.h_count))
        elif specs == (None, None, 0, 0) and (atom.element in ORGANIC_SUBSET):
            builder.append("H0")
        if atom.charge != 0:
            builder.append("{:+}".format(atom.charge))
        if brackets:
            builder.append("]")

        return "".join(builder)


def bond_to_smiles(bond: DirectedBond) -> str:
    if bond.order == 1:
        return bond.stereo if (bond.stereo in SMILES_STEREO_BONDS) else ""
    elif bond.order == 2:
        return "="
    elif bond.order == 3:
        return "#"
    else:  # this should never happen
        raise ValueError()


def mol_to_smiles(mol: MolecularDFSTree) -> str:
    assert mol.is_kekulized()

    fragments = []
    ring_log = dict()
    for root in mol.get_roots():
        derived = []
        _derive_smiles_from_fragment(derived, mol, root, ring_log)
        fragments.append("".join(derived))
    return ".".join(fragments)


def _derive_smiles_from_fragment(derived, mol, root, ring_log):
    curr_atom, curr = mol.get_atom(root), root
    derived.append(atom_to_smiles(curr_atom))

    out_bonds = mol.get_out_dirbonds(curr)
    for i, bond in enumerate(out_bonds):
        if bond.ring_bond:
            derived.append(bond_to_smiles(bond))
            ends = (min(bond.src, bond.dst), max(bond.src, bond.dst))
            rnum = ring_log.setdefault(ends, len(ring_log) + 1)
            if rnum >= 10:
                derived.append("%")
            derived.append(str(rnum))

        else:
            if i < len(out_bonds) - 1:
                derived.append("(")

            derived.append(bond_to_smiles(bond))
            _derive_smiles_from_fragment(derived, mol, bond.dst, ring_log)

            if i < len(out_bonds) - 1:
                derived.append(")")
