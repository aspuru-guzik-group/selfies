from typing import Optional

from rdkit.Chem import Kekulize, MolFromSmiles, MolToSmiles

from selfies.grammar_rules import get_chars_from_n, get_num_from_bond


def encoder(smiles: str, print_error: bool = False) -> Optional[str]:
    """Translates a SMILES into a SELFIES.

    The SMILES to SELFIES translation occurs independently of the SELFIES
    alphabet and grammar. Thus, ``selfies.encoder`` will work regardless of
    the alphabet and grammar rules that ``selfies`` is operating on, assuming
    the input is a valid SMILES. Additionally, ``selfies.encoder`` maintains the
    atom and branch order of the input SMILES; thus, one could generate random
    SELFIES corresponding to the same molecule by generating random SMILES,
    and then translating them.

    However, encoding and then decoding a SMILES may not necessarily yield
    the original SMILES. Reasons include:
        1.  SMILES with explicit aromatic symbols are automatically
            Kekulized before being translated.
        2.  A SMILES that violates the bond constraints specified by ``selfies``
            will be successfully encoded by ``selfies.encoder``, but then
            decoded into a new molecule that satisfies the constraints.
        3.  The exact ring numbering order is lost in ``selfies.encoder``,
            and cannot be reconstructed by ``selfies.decoder``.

    Finally, note that ``selfies.encoder`` does **not** check if the input
    SMILES is valid, and should not be expected to reject invalid inputs.

    :param smiles: The SMILES to be translated.
    :param print_error: If True, error messages will be printed to console.
        Defaults to False.
    :return: the SELFIES translation of **smiles**. If an error occurs,
        and **smiles** cannot be translated, :code:`None` is returned instead.

    :Example:

    >>> import selfies
    >>> selfies.encoder('C=CF')
    '[C][=C][F]'

    .. warning:: Currently, :code:`selfies.encoder` first splits a
        SMILES by the dot-bond character, translates each fragment separately,
        and then rejoins them. Thus, SMILES such as C1.C2.C12 (propane) and
        c1cc([O-].[Na+])ccc1 (sodium phenoxide) are **not** supported because
        their fragments (e.g. C1, c1cc([O-], C12) are not valid SMILES.
        However, [Na+].[O-]c1ccccc1 and [K].Cl are supported.

    """

    # TODO: remove me later
    if 'c' in smiles:
        m = MolFromSmiles(smiles)
        if m is None:
            return None

        Kekulize(m)
        smiles = MolToSmiles(m, kekulizeSmiles=True)
    # TODO: remove me later

    try:
        smiles.replace('-', '')  # remove explicit single bonds

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


def _parse_smiles(smiles):
    """A generator, which parses a SMILES string and returns character(s) of it
    one-by-one as a tuple of:
        (1) the bond character connecting the current character to the previous
            SMILES character (e.g. '=', '', '#')
        (2) the character(s) as a string (e.g. 'C', '12', '(')
        (3) the type of character(s), represented as an integer that is either
            <ATOM_TYPE>, <BRANCH_TYPE>, and <RING_TYPE>.
    As a precondition, we also assume <smiles> has no dots in it.

    Args:
        smiles: the SMILES string (without dots) to be parsed

    Returns:
        str: the character(s) of <smiles> along with their types
    """

    i = 0

    while 0 <= i < len(smiles):

        bond = ''

        if smiles[i] in {'/', '\\', '=', '#'}:
            bond = smiles[i]
            i += 1
        elif smiles[i] == '-':
            i += 1  # TODO: check if ignoring is always valid

        if 'A' <= smiles[i] <= 'Z' or smiles[i] == '*':  # elements or wildcard
            if 'a' <= smiles[i + 1: i + 2] <= 'z':  # two letter elements
                char = smiles[i: i + 2]
                char_type = ATOM_TYPE
                i += 2
            else:
                char = smiles[i]  # one letter elements (e.g. C, N, ...)
                char_type = ATOM_TYPE
                i += 1

        elif smiles[i] in ['(', ')']:  # open and closed branch brackets
            bond = smiles[i + 1]
            char = smiles[i]
            char_type = BRANCH_TYPE
            i += 1

        elif smiles[i] == '[':  # atoms encased in brackets (e.g. [NH])
            r_idx = smiles.find(']', i + 1)
            char = smiles[i + 1: r_idx] + "expl"
            char_type = ATOM_TYPE
            i = r_idx + 1

        elif '0' <= smiles[i] <= '9':  # one-digit ring number
            char = smiles[i]
            char_type = RING_TYPE
            i += 1

        elif smiles[i] == '%':  # two-digit ring number (e.g. %12)
            char = smiles[i + 1: i + 3]
            char_type = RING_TYPE
            i += 3

        else:
            raise ValueError(f"Unknown symbol '{smiles[i]}' in SMILES.")

        yield bond, char, char_type


def _translate_smiles(smiles):
    """A helper for selfies.encoder, which converts a SMILES string
    without dots into its SELFIES representation.

    Args:
        smiles: the SMILES string (without dots) to be encoded

    Returns:
        str: the SELFIES translation of <smiles>.
    """

    smiles_gen = _parse_smiles(smiles)

    # a simple mutable counter to track which atom was the i-th derived atom
    derive_counter = [0]

    # a dictionary to keep track of the rings to be made. If a ring with id
    # X is connected to the i-th and j-th derived atoms (i < j) with bond
    # character s, then after the i-th atom is derived, rings[X] = (s, i).
    # As soon as the j-th atom is derived, rings[X] is removed from <rings>,
    # and the ring is made.
    rings = {}

    selfies, _ = _translate_smiles_derive(smiles_gen, derive_counter, rings)

    return selfies


def _translate_smiles_derive(smiles_gen, counter, rings):
    """Recursive helper for _translate_smiles. Derives the SELFIES characters
    from a SMILES string, and returns a tuple of (1) the encoded SELFIES
    and (2) the length of the encoded SELFIES (as in the number of SELFIES
    characters in the string).

    Args:
        smiles_gen: a generator produced by calling _parse_smiles on a
                    SMILES string.
        counter: see <derived_counter> in _translate_smiles
        rings: see <rings> in _translate_smiles

    Returns:
        tuple: a tuple of the encoded SELFIES and its length
    """

    selfies = ""
    selfies_len = 0

    for i, (bond, char, char_type) in enumerate(smiles_gen):

        if char_type == ATOM_TYPE:
            selfies += f"[{bond}{char}]"
            counter[0] += 1
            selfies_len += 1

        elif char_type == BRANCH_TYPE:
            if char == '(':

                branch, branch_len = \
                    _translate_smiles_derive(smiles_gen, counter, rings)

                N_as_chars = get_chars_from_n(branch_len - 1)
                bond_num = get_num_from_bond(bond)

                selfies += f"[Branch{len(N_as_chars)}_{bond_num}]"
                selfies += ''.join(N_as_chars) + branch
                selfies_len += 1 + len(N_as_chars) + branch_len

            else:  # char == ')'
                return selfies, selfies_len

        else:  # char_type == RING_TYPE
            ring_id = int(char)

            if ring_id in rings:
                left_bond, left_end = rings.pop(ring_id)
                right_bond, right_end = bond, counter[0]

                ring_len = right_end - left_end
                N_as_chars = get_chars_from_n(ring_len - 1)

                if left_bond != '':
                    selfies += f"[Expl{left_bond}Ring{len(N_as_chars)}]"
                elif right_bond != '':
                    selfies += f"[Expl{right_bond}Ring{len(N_as_chars)}]"
                else:
                    selfies += f"[Ring{len(N_as_chars)}]"

                selfies += ''.join(N_as_chars)
                selfies_len += 1 + len(N_as_chars)

            else:
                rings[ring_id] = (bond, counter[0])

    return selfies, selfies_len
