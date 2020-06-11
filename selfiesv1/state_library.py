import itertools
from typing import Dict

_valences = {
    'H': 1, 'F': 1, 'Cl': 1, 'Br': 1, 'I': 1,
    'O': 2, '[NH]': 2,
    'N': 3, 'P': 3, '[C@H]': 3, '[C@@H]': 3,
    'C': 4, '[C@]': 4, '[C@@]': 4,
    'S': 6
}

_bonds = {'': 1, '/': 1, '\\': 1, '=': 2, '#': 3}


def build_state_dict() -> Dict:
    state_library = {}

    for i in range(4):

        state_dict = dict()
        state_dict['[epsilon]'] = ('', 0) if i == 0 else ('', -1)

        for (atom, valence), (bond, bond_num) in \
                itertools.product(_valences.items(), _bonds.items()):

            if bond_num > valence:
                continue

            if atom[0] == '[':
                selfies_atom = f"[{bond}{atom[1:-1]}expl]"
            else:
                selfies_atom = f"[{bond}{atom}]"

            if atom == 'H':  # edge case
                atom = '[H]'

            if i == 0:
                state_dict[selfies_atom] = (atom, valence)
                continue

            if bond_num > i:
                bond = ('', '=', '#')[i - 1]
                bond_num = i

            next_state = valence - bond_num
            if next_state == 0:
                next_state = -1

            state_dict[selfies_atom] = (bond + atom, next_state)

        state_dict['[???]'] = (None, 6)
        state_library[i] = state_dict

    state_library[4] = state_library[3]
    state_library[5] = state_library[3]
    state_library[6] = state_library[3]
    state_library[9991] = state_library[1]
    state_library[9992] = state_library[2]
    state_library[9993] = state_library[3]

    return state_library
