import itertools
from typing import Dict

# key = SMILES atom, value = maximum number of bonds key can make
_bond_orders = {
    'H': 1, 'F': 1, 'Cl': 1, 'Br': 1, 'I': 1,
    'O': 2, '[NH]': 2,
    'N': 3, '[C@H]': 3, '[C@@H]': 3,
    'C': 4, '[C@]': 4, '[C@@]': 4,
    'P': 5,
    'S': 6
}

# key = bond type, value = bond multiplicity
# '-' (explicit) single bond is not included because it is always ignored
# by the SELFIES encoder
_bonds = {'': 1, '/': 1, '\\': 1, '=': 2, '#': 3}


def build_state_dict() -> Dict:
    """Dynamically builds the state dict used to enforce the SELFIES
    grammar, based on <_bond_orders> and <_bonds>. See <_state_library>
    in utils.py for further explanation.

    Returns: the SELFIES derivation state dict
    """

    state_library = {}

    for state in range(4):

        state_dict = dict()
        state_dict['[epsilon]'] = ('', 0) if state == 0 else ('', -1)

        for (atom, bond_order), (bond, bond_num) in \
                itertools.product(_bond_orders.items(), _bonds.items()):

            if bond_num > bond_order:
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

                next_state = bond_order - bond_num
                if next_state == 0:
                    next_state = -1

                state_dict[selfies_atom] = (bond + atom, next_state)

            else:
                state_dict[selfies_atom] = (atom, bond_order)

        state_dict['[???]'] = (None, 6)
        state_library[state] = state_dict

    state_library[4] = state_library[3]
    state_library[5] = state_library[3]
    state_library[6] = state_library[3]
    state_library[9991] = state_library[1]
    state_library[9992] = state_library[2]
    state_library[9993] = state_library[3]

    return state_library
