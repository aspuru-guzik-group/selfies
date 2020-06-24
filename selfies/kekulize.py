import pprint

ATOM_TYPE = 1
BRANCH_TYPE = 2
RING_TYPE = 3


def kekulize_parser(smiles_gen):
    smiles_chars = list(map(list, smiles_gen))

    mol_graph = AromaticGraph(smiles_chars)
    _build_molecular_graph(mol_graph, smiles_chars)
    _kekulize(mol_graph)

    pprint.pprint(mol_graph.mol_graph)

    for x in mol_graph.smiles_chars:
        yield tuple(x)


def _build_molecular_graph(graph, smiles_chars, prev_idx=-1, curr_idx=-1):
    rings = {}

    while curr_idx + 1 < len(smiles_chars):

        curr_idx += 1
        _, char, char_type = smiles_chars[curr_idx]

        if char_type == ATOM_TYPE:
            if prev_idx >= 0:
                graph.add_edge(prev_idx, curr_idx, curr_idx)
            prev_idx = curr_idx

        elif char_type == BRANCH_TYPE:
            if char == '(':
                curr_idx = _build_molecular_graph(graph, smiles_chars,
                                                  prev_idx, curr_idx)
            else:
                break

        else:
            if char in rings:
                graph.add_edge(prev_idx, rings.pop(char), curr_idx)
            else:
                rings[char] = prev_idx

    return curr_idx


def _kekulize(mol_graph):
    pass


# Aromatic Helper Methods and Classes

_aromatics = {'c', 'o', 's', 'n', 'n+expl', 'n-expl'}


def _is_aromatic_char(char):
    if char.lower() == char and char not in _aromatics:
        raise ValueError(f"Kekulization of '{char}' is not supported.")

    return char in _aromatics


class AromaticGraph:

    def __init__(self, smiles_chars):
        self.smiles_chars = smiles_chars
        self.mol_graph = {}

    def get_atom_char(self, idx):
        return self.smiles_chars[idx][1]

    def get_bond_char(self, idx):
        return self.smiles_chars[idx][0]

    def get_directed_edge(self, idx_a, idx_b):

        if idx_a not in self.mol_graph:
            return None

        for edge in self.mol_graph[idx_a]:
            if edge[0] == idx_b:
                return edge

        return None

    def add_edge(self, idx_a, idx_b, bond_idx):

        is_a_aromatic = _is_aromatic_char(self.get_atom_char(idx_a))
        is_b_aromatic = _is_aromatic_char(self.get_atom_char(idx_b))
        bond_char = self.get_bond_char(bond_idx)

        if is_a_aromatic or is_b_aromatic:

            if is_a_aromatic and is_b_aromatic and bond_char == '':
                bond_char = ':'  # aromatic bond

            if is_a_aromatic:
                a_edges = self.mol_graph.setdefault(idx_a, [])
                a_edges.append([idx_b, bond_char, bond_idx])

            if is_b_aromatic:
                b_edges = self.mol_graph.setdefault(idx_b, [])
                b_edges.append([idx_a, bond_char, bond_idx])

    def update_edge(self, idx_a, idx_b, new_bond_char):

        a_to_b = self.get_directed_edge(idx_a, idx_b)
        b_to_a = self.get_directed_edge(idx_b, idx_a)

        a_to_b[1] = new_bond_char
        b_to_a[1] = new_bond_char


if __name__ == '__main__':
    from selfies.encoder import _parse_smiles

    list(kekulize_parser(_parse_smiles("[H]n1cccc1")))
