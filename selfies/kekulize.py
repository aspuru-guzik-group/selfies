import pprint

ATOM_TYPE = 1
BRANCH_TYPE = 2
RING_TYPE = 3


def kekulize_parser(smiles_gen):
    smiles_chars = list(map(list, smiles_gen))

    mol_graph = AromaticGraph(smiles_chars)

    rings = {}
    _build_molecular_graph(mol_graph, smiles_chars, rings)

    _kekulize(mol_graph)

    m = mol_graph.mol_graph
    for c in m:
        m[c] = list(map(str, m[c]))
    pprint.pprint(m)

    for x in mol_graph.smiles_chars:
        yield tuple(x)


def _build_molecular_graph(graph, smiles_chars, rings,
                           prev_idx=-1, curr_idx=-1):
    while curr_idx + 1 < len(smiles_chars):

        curr_idx += 1
        _, char, char_type = smiles_chars[curr_idx]

        if char_type == ATOM_TYPE:
            if prev_idx >= 0:
                graph.add_edge(prev_idx, curr_idx, curr_idx)
            prev_idx = curr_idx

        elif char_type == BRANCH_TYPE:
            if char == '(':
                curr_idx = _build_molecular_graph(graph, smiles_chars, rings,
                                                  prev_idx, curr_idx)
            else:
                break

        else:
            if char in rings:
                left_idx, left_bond_idx = rings.pop(char)
                right_idx, right_bond_idx = prev_idx, curr_idx

                if smiles_chars[left_bond_idx][0] != '':
                    bond_idx = left_bond_idx
                    smiles_chars[right_bond_idx][0] = ''
                else:
                    bond_idx = right_bond_idx
                    smiles_chars[left_bond_idx][0] = ''

                graph.add_edge(left_idx, right_idx, bond_idx)
            else:
                rings[char] = (prev_idx, curr_idx)

    return curr_idx


def _kekulize(mol_graph):
    mol_graph.prune()


# Aromatic Helper Methods and Classes

_aromatics = {
    'c': 4,
    'o': 6,
    's': 6,
    'n': 5,
    'n+expl': 4,
    'n-expl': 6,
}

_double_bonds = {
    3: [1, 0],
    4: [1, 1],
    5: [1, 0],
    6: [0, 1]
}


def _is_aromatic_char(char):
    if char.lower() == char and char not in _aromatics:
        raise ValueError(f"Kekulization of '{char}' is not supported.")

    return char in _aromatics


class AromaticGraph:

    def __init__(self, smiles_chars):
        self.smiles_chars = smiles_chars
        self.mol_graph = {}
        self.aromatic_atoms = set()

    def get_atom_char(self, idx):
        return self.smiles_chars[idx][1]

    def get_bond_char(self, idx):
        return self.smiles_chars[idx][0]

    def add_edge(self, idx_a, idx_b, bond_idx):

        atom_a = self.get_atom_char(idx_a)
        atom_b = self.get_atom_char(idx_b)
        bond_char = self.get_bond_char(bond_idx)

        if _is_aromatic_char(atom_a):
            self.aromatic_atoms.add(idx_a)

        if _is_aromatic_char(atom_b):
            self.aromatic_atoms.add(idx_b)

        if bond_char == ':':
            self.aromatic_atoms.add(idx_a)
            self.aromatic_atoms.add(idx_b)
            bond_char = ''

        edge = Edge(idx_a, idx_b, bond_char, bond_idx)

        self.mol_graph.setdefault(idx_a, []).append(edge)
        self.mol_graph.setdefault(idx_b, []).append(edge)

    def prune(self):
        non_aromatic = self.mol_graph.keys() - self.aromatic_atoms
        for idx in non_aromatic:
            self.mol_graph.pop(idx)

        print(non_aromatic)


class Edge:

    def __init__(self, idx_a, idx_b, bond_char, bond_idx):
        self.idx_a = idx_a
        self.idx_b = idx_b
        self.bond_char = bond_char
        self.bond_idx = bond_idx

    def __str__(self):  # for debugging
        bond = self.bond_char
        if bond == '':
            bond = '-'

        return f"[{self.idx_a}{bond}{self.idx_b} | {self.bond_idx}]"


if __name__ == '__main__':
    from selfies.encoder import _parse_smiles

    list(kekulize_parser(_parse_smiles("CC[H]n1cccc=1")))
