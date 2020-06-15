import itertools
from typing import Dict

# default atom_dict
# key = SMILES atom, value = maximum number of bonds key can make
default_atom_dict = {
    'H': 1, 'F': 1, 'Cl': 1, 'Br': 1, 'I': 1,
    'O': 2, '[NH]': 2,
    'N': 3, '[C@H]': 3, '[C@@H]': 3,
    'C': 4, '[C@]': 4, '[C@@]': 4,
    'P': 5,
    'S': 6,
}

# key = bond type, value = bond multiplicity
# '-' (explicit) single bond is not included because it is always ignored
# by the SELFIES encoder
_bonds = {'': 1, '/': 1, '\\': 1, '=': 2, '#': 3}


def build_state_dict(atom_dict=None) -> Dict:
    """Dynamically builds the state dict used to enforce the SELFIES
    grammar based on <atom_dict>. See <state_library> in utils.py for
    further explanation of the structure of the state dict.

    Args:
        atom_dict: a dictionary of the maximum bond capacity of various atom(s)
                   or ions, as specified by <set_selfies_alphabet> in
                   utils.py. If None, then has value <default_atom_dict>.

    Returns: the SELFIES derivation state dict
    """

    if atom_dict is None:
        atom_dict = default_atom_dict

    state_library = {}

    for state in range(4):

        state_dict = dict()
        state_dict['[epsilon]'] = ('', -1)

        for (atom, max_bonds), (bond, bond_num) in \
                itertools.product(atom_dict.items(), _bonds.items()):

            if bond_num > max_bonds:
                continue  # e.g. can't make a double bond with F

            if atom[0] == '[':
                selfies_atom = f"[{bond}{atom[1:-1]}expl]"
            else:
                selfies_atom = f"[{bond}{atom}]"

            if atom == 'H':  # edge case: [H] (SELFIES) --> [H] (SMILES)
                atom = '[H]'

            if state > 0:

                if bond_num > state:
                    bond = ('', '=', '#')[state - 1]
                    bond_num = state

                next_state = max_bonds - bond_num
                if next_state == 0:
                    next_state = -1

                state_dict[selfies_atom] = (bond + atom, next_state)

            else:
                state_dict[selfies_atom] = (atom, max_bonds)

        state_dict['[???]'] = (None, 8)
        state_library[state] = state_dict

    state_library[4] = state_library[3]
    state_library[5] = state_library[3]
    state_library[6] = state_library[3]
    state_library[7] = state_library[3]
    state_library[8] = state_library[3]
    state_library[9991] = state_library[1]
    state_library[9992] = state_library[2]
    state_library[9993] = state_library[3]

    return state_library
