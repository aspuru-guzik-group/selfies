from selfiesv1.state_dicts import get_bond_num, index_to_chars


def encoder(smiles, print_error=True):
    try:
        smiles = smiles.replace(" ", '')

        all_smiles = []  # process dot-separated fragments separately
        for s in smiles.split("."):
            all_smiles.append(_translate_smiles(s))
        return '.'.join(all_smiles)

    except ValueError as err:
        if print_error:
            print(err)
            print('Could not encode SMILES. Please contact authors.')
        return -1


ATOM_TYPE = 1
BRANCH_TYPE = 2
RING_TYPE = 3


def _parse_smiles(smiles):
    i = 0

    while 0 <= i < len(smiles):

        bond = ''

        if smiles[i] in {'/', '\\', '=', '#'}:
            bond = smiles[i]
            i += 1
        elif smiles[i] == '-':
            i += 1  # TODO: check if ignoring is always valid

        if 'A' <= smiles[i] <= 'Z' or smiles[i] == '*':
            if smiles[i: i + 2] in {'Br', 'Cl'}:
                char = (smiles[i: i + 2], ATOM_TYPE)
                i += 2
            else:
                char = (smiles[i], ATOM_TYPE)
                i += 1

        elif smiles[i] in ['(', ')']:
            bond = smiles[i + 1]
            char = (smiles[i], BRANCH_TYPE)
            i += 1

        elif smiles[i] == '[':
            r_idx = smiles.find(']', i + 1)
            char = (smiles[i + 1: r_idx] + "expl", ATOM_TYPE)
            i = r_idx + 1

        elif '0' <= smiles[i] <= '9':
            char = (smiles[i], RING_TYPE)
            i += 1

        elif smiles[i] == '%':
            char = (smiles[i + 1: i + 3], RING_TYPE)
            i += 3

        else:
            raise ValueError("Unknown symbol in SMILES.")

        yield bond, char


def _translate_smiles(smiles):
    smiles_gen = _parse_smiles(smiles)
    derive_counter = [0]
    rings = {}

    selfies, _ = _translate_smiles_derive(smiles_gen, derive_counter, rings)

    return selfies


def _translate_smiles_derive(smiles_gen, counter, rings):
    selfies = ""
    selfies_len = 0

    for i, (bond, (char, char_type)) in enumerate(smiles_gen):

        if char_type == ATOM_TYPE:
            selfies += f"[{bond}{char}]"
            counter[0] += 1
            selfies_len += 1

        elif char_type == BRANCH_TYPE:
            if char == '(':

                branch, branch_len = \
                    _translate_smiles_derive(smiles_gen, counter, rings)

                N_as_chars = index_to_chars(branch_len - 1)
                bond_num = get_bond_num(bond)

                selfies += f"[Branch{len(N_as_chars)}_{bond_num}]"
                selfies += ''.join(N_as_chars) + branch
                selfies_len += 1 + len(N_as_chars) + branch_len

            else:
                return selfies, selfies_len

        else:
            ring_id = int(char)

            if ring_id in rings:
                left_bond, left_end = rings.pop(ring_id)
                right_bond, right_end = bond, counter[0]

                ring_len = right_end - left_end
                N_as_chars = index_to_chars(ring_len - 1)

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
