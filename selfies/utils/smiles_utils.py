import enum
import re
from collections import deque
from typing import Iterator, Optional, Tuple, Union, List

from selfies.constants import AROMATIC_SUBSET, ELEMENTS, ORGANIC_SUBSET
from selfies.exceptions import SMILESParserError
from selfies.mol_graph import Atom, Attribution, \
    AttributionMap, DirectedBond, MolecularGraph

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
    """A token in a SMILES string, containing a symbol (atom, branch bracket,
    ring number, dot) and its preceding bond, if it exists (e.g. =C, %12, #N).
    """

    def __init__(
            self,
            bond_idx: Optional[int],
            start_idx: int, end_idx: int, token_type: SMILESTokenTypes,
            token: str
    ):
        self.bond_idx = bond_idx
        self.start_idx = start_idx
        self.end_idx = end_idx
        self.token_type = token_type
        self.token = token

    def extract_bond_char(self, smiles):
        return None if (self.bond_idx is None) else smiles[self.bond_idx]

    def extract_symbol(self, smiles):
        return smiles[self.start_idx:self.end_idx]

    def __str__(self):
        return self.token


def tokenize_smiles(smiles: str) -> Iterator[SMILESToken]:
    """Splits a SMILES string into its tokens.

    :param smiles: the input SMILES string.
    :return: the tokens of the input SMILES one-by-one with order preserved.
    """

    i = 0
    while i < len(smiles):

        if smiles[i] == ".":
            yield SMILESToken(None, i, i + 1, SMILESTokenTypes.DOT, smiles[i])
            i += 1
            continue

        if smiles[i] in SMILES_BOND_ORDERS:
            bond_idx = i
            i += 1
        else:
            bond_idx = None

        if i == len(smiles):
            raise SMILESParserError(smiles, "hanging bond", i - 1)

        elif smiles[i].isalpha():  # organic subset elements
            if smiles[i: i + 2] in ("Br", "Cl"):  # two-letter elements
                token = SMILESToken(bond_idx, i, i + 2,
                                    SMILESTokenTypes.ATOM, smiles[i: i + 2])
            else:  # one-letter elements (e.g. C, N, ...)
                token = SMILESToken(bond_idx, i, i + 1,
                                    SMILESTokenTypes.ATOM, smiles[i:i + 1])

        elif smiles[i] == "[":  # atoms encased in brackets (e.g. [NH])
            r_idx = smiles.find("]", i + 1)
            if r_idx == -1:
                raise SMILESParserError(smiles, "hanging bracket [", i)
            token = SMILESToken(bond_idx, i, r_idx + 1,
                                SMILESTokenTypes.ATOM, smiles[i:r_idx + 1])

        elif smiles[i] in ("(", ")"):  # open and closed branch brackets
            if bond_idx is not None:
                raise SMILESParserError(smiles, "hanging_bond", bond_idx)
            token = SMILESToken(
                None, i, i + 1, SMILESTokenTypes.BRANCH, smiles[i:i+1])

        elif smiles[i].isdigit():  # one-digit ring number
            token = SMILESToken(bond_idx, i, i + 1,
                                SMILESTokenTypes.RING, smiles[i:i+1])

        elif smiles[i] == "%":  # two-digit ring number (e.g. %12)
            rnum = smiles[i + 1: i + 3]
            if not (rnum.isnumeric() and len(rnum) == 2):
                err_msg = "invalid ring number '%{}'".format(rnum)
                raise SMILESParserError(smiles, err_msg, i)
            token = SMILESToken(bond_idx, i, i + 3,
                                SMILESTokenTypes.RING, smiles[i:i+3])

        else:
            err_msg = "unrecognized symbol '{}'".format(smiles[i])
            raise SMILESParserError(smiles, err_msg, i)

        yield token
        i = token.end_idx


# =============================================================================
# SMILES -> Atom, Graph, etc.
# =============================================================================


def smiles_to_atom(atom_symbol: str) -> Optional[Atom]:
    """Reads an atom from its SMILES representation.

    :param atom_symbol: a SMILES atom symbol.
    :return: the atom that the input symbol represents.
    """

    if atom_symbol[0] == "[" and atom_symbol[-1] == "]":
        pass  # continue below
    elif atom_symbol in ORGANIC_SUBSET:  # e.g. C, N, O, ...
        return Atom(atom_symbol, False)
    elif atom_symbol in AROMATIC_SUBSET:  # e.g. c, n, o, ...
        return Atom(atom_symbol.capitalize(), True)
    else:
        return None

    # e.g. [C], [C@@H], [O-], ...
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
    """Reads a bond from its SMILES representation.

    :param bond_char: a SMILES bond symbol.
    :return: the order and stereochemical specification of the bond
        that the input symbol represents.
    """

    order = SMILES_BOND_ORDERS.get(bond_char, 1)
    stereo = bond_char if (bond_char in SMILES_STEREO_BONDS) else None
    return order, stereo


def smiles_to_mol(smiles: str, attributable: bool) -> MolecularGraph:
    """Reads a molecular graph from a SMILES string.

    :param smiles: the input SMILES string.
    :param attributable: if molecular graph needs to include attributions
    :return: a molecular graph that the input SMILES string represents.
    :raises SMILESParserError: if the input SMILES is invalid.
    """

    if smiles == "":
        raise SMILESParserError(smiles, "empty SMILES", 0)

    mol = MolecularGraph(attributable=attributable)
    tokens = deque(tokenize_smiles(smiles))
    i = 0
    while tokens:
        i = _derive_mol_from_tokens(mol, smiles, tokens, i)
    return mol


def _derive_mol_from_tokens(mol, smiles, tokens, i):
    tok = None
    prev_stack = deque()  # keep track of previous atom on the current chain
    branch_stack = deque()  # keep track of open branches
    ring_log = dict()  # keep track of hanging ring numbers
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

            curr, i = _attach_atom(mol, bond_char, curr, prev_atom, i, tok)
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
        i += 1

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
    return i


def _attach_atom(mol, bond_char, atom, prev_atom, i, tok):
    is_root = (prev_atom is None)
    if bond_char:
        i += 1
    o = mol.add_atom(atom, mark_root=is_root)
    mol.add_attribution(o, [Attribution(i, str(tok))])
    if not is_root:
        src, dst = prev_atom.index, atom.index
        order, stereo = smiles_to_bond(bond_char)
        if prev_atom.is_aromatic and atom.is_aromatic and (bond_char is None):
            order = 1.5  # handle implicit aromatic bonds, e.g. cc
        o = mol.add_bond(src=src, dst=dst, order=order, stereo=stereo)
        mol.add_attribution(o, [Attribution(i, str(tok))])
    return atom, i


def _make_ring_bonds(mol, smiles, ltoken, latom, lpos, rtoken, ratom):
    if mol.has_bond(latom.index, ratom.index):
        err_msg = "ring bond specified between already-bonded atoms"
        raise SMILESParserError(smiles, err_msg, ltoken.start_idx)

    lbond_char = ltoken.extract_bond_char(smiles)
    rbond_char = rtoken.extract_bond_char(smiles)

    # checking that ring bonds match
    bonds = (lbond_char, rbond_char)
    if bonds[0] is None:
        bonds = (bonds[1], bonds[0])
    # swap bonds so that if bonds[0] is None, then bonds[1] is None

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
        lorder = rorder = 1.5  # handle implicit aromatic bonds, e.g. c1ccccc1

    mol.add_ring_bond(
        a=latom.index, a_stereo=lstereo, a_pos=lpos,
        b=ratom.index, b_stereo=rstereo,
        order=max(lorder, rorder)
    )


# =============================================================================
# SMILES <- Atom, Graph, etc.
# =============================================================================


def atom_to_smiles(atom: Atom, brackets: bool = True) -> str:
    """Converts an atom into its SMILES representation.

    :param atom: the input atom.
    :param brackets: True, if brackets should be added around the returned
        symbol (e.g. in the case of [C] or [C@@H]). Defaults to True.
    :return: a SMILES symbol representing the input atom.
    """
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
    """Converts a bond into its SMILES representation.

    :param bond: the input bond.
    :return: a SMILES symbol representing the input bond.
    """

    if bond.order == 1:
        return bond.stereo if (bond.stereo in SMILES_STEREO_BONDS) else ""
    elif bond.order == 2:
        return "="
    elif bond.order == 3:
        return "#"
    else:  # this should never happen
        raise ValueError()


def mol_to_smiles(
    mol: MolecularGraph,
    attribute: bool = False
) -> Union[str, Tuple[str, List[Tuple[str,  List[Tuple[int, str]]]]]]:
    """Converts a molecular graph into its SMILES representation, maintaining
    the traversal order indicated by the input graph.

    :param mol: the input molecule.
    :param attribute: if an attribution should be returned
    :return: a SMILES string representing the input molecule if
             attribute is ``False``, otherwise a tuple is returned of
             SMILES string and attribution list.
    """
    assert mol.is_kekulized()

    fragments = []
    attribution_maps = []
    attribution_index = 0
    ring_log = dict()
    for root in mol.get_roots():
        derived = []
        _derive_smiles_from_fragment(
            derived, mol, root, ring_log, attribution_maps, attribution_index)
        attribution_index += _strlen(derived)
        fragments.append("".join(derived))
    # trim attribution map of empty tokens
    attribution_maps = [a for a in attribution_maps if a.token]
    result = ".".join(fragments), attribution_maps
    return result if attribute else result[0]


def _strlen(slist: List[str]) -> int:
    return len(''.join(slist))


def _derive_smiles_from_fragment(
        derived,
        mol,
        root,
        ring_log,
        attribution_maps, attribution_index=0):
    curr_atom, curr = mol.get_atom(root), root
    token = atom_to_smiles(curr_atom)
    derived.append(token)
    attribution_maps.append(AttributionMap(
        _strlen(derived) - 1 + attribution_index,
        token, mol.get_attribution(curr_atom)))

    out_bonds = mol.get_out_dirbonds(curr)
    for i, bond in enumerate(out_bonds):
        if bond.ring_bond:
            token = bond_to_smiles(bond)
            derived.append(token)
            attribution_maps.append(AttributionMap(
                _strlen(derived) - 1 + attribution_index,
                token, mol.get_attribution(bond)))
            ends = (min(bond.src, bond.dst), max(bond.src, bond.dst))
            rnum = ring_log.setdefault(ends, len(ring_log) + 1)
            if rnum >= 10:
                derived.append("%")
            derived.append(str(rnum))

        else:
            if i < len(out_bonds) - 1:
                derived.append("(")

            token = bond_to_smiles(bond)
            derived.append(token)
            attribution_maps.append(AttributionMap(
                _strlen(derived) - 1 + attribution_index,
                token, mol.get_attribution(bond)))
            _derive_smiles_from_fragment(
                derived, mol, bond.dst, ring_log,
                attribution_maps, attribution_index)
            if i < len(out_bonds) - 1:
                derived.append(")")
    return attribution_maps
