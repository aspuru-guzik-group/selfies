from typing import Dict, Iterable, List, Optional, Tuple

from selfies.grammar_rules import get_num_from_bond, get_symbols_from_n
from selfies.kekulize import kekulize_parser


def encoder(smiles: str, print_error: bool = False) -> Optional[str]:
    """Translates a SMILES into a SELFIES.

    The SMILES to SELFIES translation occurs independently of the SELFIES
    alphabet and grammar. Thus, :func:`selfies.encoder` will work regardless of
    the alphabet and grammar rules that :py:mod:`selfies` is operating on,
    assuming the input is a valid SMILES. Additionally, :func:`selfies.encoder`
    preserves the atom and branch order of the input SMILES; thus, one
    could generate random SELFIES corresponding to the same molecule by
    generating random SMILES, and then translating them.

    However, encoding and then decoding a SMILES may not necessarily yield
    the original SMILES. Reasons include:

        1.  SMILES with aromatic symbols are automatically
            Kekulized before being translated.
        2.  SMILES that violate the bond constraints specified by
            :mod:`selfies` will be successfully encoded by
            :func:`selfies.encoder`, but then decoded into a new molecule
            that satisfies the constraints.
        3.  The exact ring numbering order is lost in :func:`selfies.encoder`,
            and cannot be reconstructed by :func:`selfies.decoder`.

    Finally, note that :func:`selfies.encoder` does **not** check if the input
    SMILES is valid, and should not be expected to reject invalid inputs.
    It is recommended to use RDKit to first verify that the SMILES are
    valid.

    :param smiles: The SMILES to be translated.
    :param print_error: If True, error messages will be printed to console.
        Defaults to False.
    :return: the SELFIES translation of ``smiles``. If an error occurs,
        and ``smiles`` cannot be translated, :code:`None` is returned instead.

    :Example:

    >>> import selfies
    >>> selfies.encoder('C=CF')
    '[C][=C][F]'

    .. note:: Currently, :func:`selfies.encoder` does not support the
        following types of SMILES:

        *   SMILES using ring numbering across a dot-bond symbol
            to specify bonds, e.g. ``C1.C2.C12`` (propane) or
            ``c1cc([O-].[Na+])ccc1`` (sodium phenoxide).
        *   SMILES with ring numbering between atoms that are over
            ``16 ** 3 = 4096`` atoms apart.
        *   SMILES using the wildcard symbol ``*``.
        *   SMILES using chiral specifications other than ``@`` and ``@@``.

    """

    try:
        all_selfies = []  # process dot-separated fragments separately
        for s in smiles.split("."):
            all_selfies.append(_translate_smiles(s))
        return '.'.join(all_selfies)

    except ValueError as err:
        if print_error:
            print(err)
            print('Could not encode SMILES. Please contact authors.')
        return None


ATOM_TYPE = 1
BRANCH_TYPE = 2
RING_TYPE = 3


def _parse_smiles(smiles: str) -> Iterable[Tuple[str, str, int]]:
    """Parses a SMILES into its symbols.

    A generator, which parses a SMILES string and returns its symbol(s)
    one-by-one as a tuple of:
        (1) the bond symbol connecting the current atom/ring/branch symbol
            to the previous atom/ring/branch symbol (e.g. '=', '', '#')
        (2) the atom/ring/branch symbol as a string (e.g. 'C', '12', '(')
        (3) the type of the symbol in (2), represented as an integer that is
            either ``ATOM_TYPE``, ``BRANCH_TYPE``, and ``RING_TYPE``.
    As a precondition, we also assume ``smiles`` has no dots in it.

    :param smiles: the SMILES to be parsed.
    :return: an iterable of the symbol(s) of the SELFIES along with
        their types.
    """

    i = 0

    while 0 <= i < len(smiles):

        bond = ''

        if smiles[i] in {'-', '/', '\\', '=', '#', ":"}:
            bond = smiles[i]
            i += 1

        if smiles[i].isalpha():  # organic subset elements
            if smiles[i: i + 2] in ('Br', 'Cl'):  # two letter elements
                symbol = smiles[i: i + 2]
                symbol_type = ATOM_TYPE
                i += 2
            else:
                symbol = smiles[i]  # one letter elements (e.g. C, N, ...)
                symbol_type = ATOM_TYPE
                i += 1

        elif smiles[i] in ('(', ')'):  # open and closed branch brackets
            bond = smiles[i + 1: i + 2]
            symbol = smiles[i]
            symbol_type = BRANCH_TYPE
            i += 1

        elif smiles[i] == '[':  # atoms encased in brackets (e.g. [NH])
            r_idx = smiles.find(']', i + 1)
            symbol = smiles[i: r_idx + 1]
            symbol_type = ATOM_TYPE
            i = r_idx + 1

        elif smiles[i].isdigit():  # one-digit ring number
            symbol = smiles[i]
            symbol_type = RING_TYPE
            i += 1

        elif smiles[i] == '%':  # two-digit ring number (e.g. %12)
            symbol = smiles[i + 1: i + 3]
            symbol_type = RING_TYPE
            i += 3

        else:
            raise ValueError("Unknown symbol '{}'.".format(smiles[i]))

        yield bond, symbol, symbol_type


def _translate_smiles(smiles: str) -> str:
    """A helper for ``selfies.encoder``, which translates a SMILES into a
    SELFIES (assuming the input SMILES contains no dots).

    :param smiles: the SMILES to be translated.
    :return: the SELFIES translation of SMILES.
    """

    smiles_gen = _parse_smiles(smiles)

    char_set = set(smiles)
    if any(c in char_set for c in ['c', 'n', 'o', 'p', 'a', 's']):
        smiles_gen = kekulize_parser(smiles_gen)

    # a simple mutable counter to track which atom was the i-th derived atom
    derive_counter = [0]

    # a dictionary to keep track of the rings to be made. If a ring with id
    # X is connected to the i-th and j-th derived atoms (i < j) with bond
    # symbol s, then after the i-th atom is derived, rings[X] = (s, i).
    # As soon as the j-th atom is derived, rings[X] is removed from <rings>,
    # and the ring is made.
    rings = {}

    selfies, _ = _translate_smiles_derive(smiles_gen, rings, derive_counter)

    return selfies


def _translate_smiles_derive(smiles_gen: Iterable[Tuple[str, str, int]],
                             rings: Dict[int, Tuple[str, int]],
                             counter: List[int]) -> Tuple[str, int]:
    """Recursive helper for _translate_smiles.

    Derives the SELFIES from a SMILES, and returns a tuple of (1) the
    translated SELFIES and (2) the symbol length of the translated SELFIES.

    :param smiles_gen: an iterable of the symbols (and their types)
        of the SMILES to be translated, created by ``_parse_smiles``.
    :param rings: See ``rings`` in ``_translate_smiles``.
    :param counter: a one-element list that serves as a mutable counter.
        See ``derived_counter`` in ``_translate_smiles``.
    :return: A tuple of the translated SELFIES and its symbol length.
    """

    selfies = ""
    selfies_len = 0
    prev_idx = -1

    for bond, symbol, symbol_type in smiles_gen:

        if bond == '-':  # ignore explicit single bonds
            bond = ''

        if symbol_type == ATOM_TYPE:
            if symbol[0] == '[':
                selfies += "[{}{}expl]".format(bond, symbol[1:-1])
            else:
                selfies += "[{}{}]".format(bond, symbol)
            counter[0] += 1
            selfies_len += 1
            prev_idx = counter[0]

        elif symbol_type == BRANCH_TYPE:
            if symbol == '(':

                # NOTE: looping inside a loop on a generator will produce
                # expected behaviour in this case.

                branch, branch_len = \
                    _translate_smiles_derive(smiles_gen, rings, counter)

                N_as_symbols = get_symbols_from_n(branch_len - 1)
                bond_num = get_num_from_bond(bond)

                selfies += "[Branch{}_{}]".format(len(N_as_symbols), bond_num)
                selfies += ''.join(N_as_symbols) + branch
                selfies_len += 1 + len(N_as_symbols) + branch_len

            else:  # symbol == ')'
                break

        else:  # symbol_type == RING_TYPE
            ring_id = int(symbol)

            if ring_id in rings:
                left_bond, left_end = rings.pop(ring_id)
                right_bond, right_end = bond, prev_idx

                ring_len = right_end - left_end
                N_as_symbols = get_symbols_from_n(ring_len - 1)

                if left_bond != '':
                    selfies += "[Expl{}Ring{}]".format(left_bond,
                                                       len(N_as_symbols))
                elif right_bond != '':
                    selfies += "[Expl{}Ring{}]".format(right_bond,
                                                       len(N_as_symbols))
                else:
                    selfies += "[Ring{}]".format(len(N_as_symbols))

                selfies += ''.join(N_as_symbols)
                selfies_len += 1 + len(N_as_symbols)

            else:
                rings[ring_id] = (bond, counter[0])

    return selfies, selfies_len
