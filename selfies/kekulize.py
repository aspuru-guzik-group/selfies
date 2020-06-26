from selfies.grammar_rules import get_num_from_bond

ATOM_TYPE = 1
BRANCH_TYPE = 2
RING_TYPE = 3


def kekulize_parser(smiles_gen):
    smiles_chars = list(map(list, smiles_gen))

    mol_graph = MolecularGraph(smiles_chars)

    rings = {}
    _build_molecular_graph(mol_graph, smiles_chars, rings)

    _kekulize(mol_graph)

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
    mol_graph.prune_to_pi_subgraph()

    visited = set()

    for i in mol_graph.graph.keys():
        matched = set()

        success = mol_graph.dfs_assign_bonds(i, visited, matched)
        if not success:
            raise ValueError("Kekulization Failed.")

    mol_graph.write_to_smiles_chars()


# Aromatic Helper Methods and Classes

_aromatic_valences = {  # wild card not supported currently
    'b': 3, 'c': 4, 'n': 5, 'p': 5, '[as]': 5,
    'o': 6, 's': 6, '[se]': 6
}


def _parse_bracketed_char(char):
    i = 1

    number = ''
    while char[i].isdigit():
        number += char[i]
        i += 1

    if char[i + 1].isalpha() and char[i + 1] != 'H':
        atom = char[i: i + 2]
        i += 2
    else:
        atom = char[i]
        i += 1

    chiral = ''
    if char[i] == '@':
        chiral = '@'
        if char[i + 1] == '@':
            chiral = '@@'
            i += 1
        i += 1

    h_count = 0
    if char[i] == 'H':
        h_count = 1

        i += 1
        if char[i].isdigit():
            h_count = int(char[i])
            i += 1

    charge = 0
    if char[i] in ('+', '-'):
        charge = 1 if char[i] == '+' else -1

        i += 1
        magnitude = ''
        while char[i].isdigit():
            magnitude += char[i]
            i += 1
        if magnitude != '':
            charge *= int(magnitude)

    class_note = char[i: -1]

    return [number, atom, chiral, h_count, charge, class_note]


def _rejoin_bracketed_char(number, atom, chiral, h_count, charge, class_note):
    # H count to character
    if h_count == 0:
        h_char = ''
    elif h_count == 1:
        h_char = 'H'
    else:
        h_char = f"H{h_count}"

    # charge to character
    if charge == 0:
        charge_char = ''
    elif abs(charge) == 1:
        charge_char = '+' if charge > 0 else '-'
    else:
        sign = '+' if charge > 0 else '-'
        charge_char = sign + abs(charge)

    return f"[{number}{atom}{chiral}{h_char}{charge_char}{class_note}]"


def _capitalize(char):
    if char[0] == '[':
        char_parts = _parse_bracketed_char(char)
        char_parts[1] = char_parts[1].capitalize()
        return _rejoin_bracketed_char(*char_parts)

    return char.upper()


def _is_aromatic(char):
    atom = char
    if char[0] == '[':
        char_parts = _parse_bracketed_char(char)
        atom = char_parts[1]

    if atom[0].isupper():
        return False

    if atom not in _aromatic_valences:
        raise ValueError(f"Kekulization Failed: aromatic symbol {char} "
                         f"not recognized.")

    return True


def _in_pi_subgraph(char, bonds):
    used_electrons = 0
    for b in bonds:
        used_electrons += get_num_from_bond(b)

    atom = char
    h_count = 0
    charge = 0

    if char[0] == '[':
        _, atom, _, h_count, charge, _ = _parse_bracketed_char(char)

    if char == 'c' and len(bonds) == 2:
        h_count += 1  # implied bonded hydrogen

    if h_count > 1:
        raise ValueError(f"Kekulization Failed: {char} not supported.")

    elif h_count == 1:  # e.g. [nH]
        used_electrons += 1

    valence = _aromatic_valences[atom] - charge
    free_electrons = valence - used_electrons
    return free_electrons % 2 != 0


class MolecularGraph:

    def __init__(self, smiles_chars):
        self.smiles_chars = smiles_chars
        self.graph = {}
        self.aro_indices = set()

    def get_atom_char(self, idx):
        return self.smiles_chars[idx][1]

    def get_bond_char(self, idx):
        return self.smiles_chars[idx][0]

    def set_bond_char(self, bond_char, idx):
        self.smiles_chars[idx][0] = bond_char

    def set_atom_char(self, atom_char, idx):
        self.smiles_chars[idx][1] = atom_char

    def add_edge(self, idx_a, idx_b, bond_idx):

        atom_a = self.get_atom_char(idx_a)
        atom_b = self.get_atom_char(idx_b)
        atom_a_aro = (idx_a in self.aro_indices) or _is_aromatic(atom_a)
        atom_b_aro = (idx_b in self.aro_indices) or _is_aromatic(atom_b)
        bond_char = self.get_bond_char(bond_idx)

        if atom_a_aro:
            self.aro_indices.add(idx_a)

        if atom_b_aro:
            self.aro_indices.add(idx_b)

        if bond_char == ':':
            self.aro_indices.add(idx_a)
            self.aro_indices.add(idx_b)

            self.set_bond_char('', bond_idx)
            bond_char = ''

        edge = Edge(idx_a, idx_b, bond_char, bond_idx)

        self.graph.setdefault(idx_a, []).append(edge)
        self.graph.setdefault(idx_b, []).append(edge)

    def prune_to_pi_subgraph(self):

        # remove nodes
        non_aromatic = self.graph.keys() - self.aro_indices
        for i in non_aromatic:
            self.graph.pop(i)

        for i in self.aro_indices:

            atom = self.get_atom_char(i)
            bonds = tuple(edge.bond_char for edge in self.graph[i])

            if not _in_pi_subgraph(atom, bonds):
                self.graph.pop(i)

        # remove irrelevant edges
        for idx, edges in self.graph.items():

            keep = list(filter(lambda e: (e.idx_a in self.graph) and
                                         (e.idx_b in self.graph) and
                                         (e.bond_char == ''),
                               edges))
            self.graph[idx] = keep

    def dfs_assign_bonds(self, idx, visited, matched):

        if idx in visited:
            return True

        edges = self.graph[idx]

        if idx in matched:
            visited.add(idx)

            for e in edges:
                adj = e.other_end(idx)
                if not self.dfs_assign_bonds(adj, visited, matched):
                    visited.remove(idx)
                    return False

            return True

        else:

            candidates = list(filter(lambda i: i.other_end(idx) not in matched,
                                     edges))

            if not candidates:
                return False

            for c in candidates:

                c.bond_char = '='
                matched.add(c.idx_a)
                matched.add(c.idx_b)

                success = self.dfs_assign_bonds(idx, visited, matched)
                if success:
                    return True
                else:
                    c.bond_char = ''
                    matched.remove(c.idx_a)
                    matched.remove(c.idx_b)

            return False

    def write_to_smiles_chars(self):

        # capitalize aromatic molecules
        for idx in self.aro_indices:
            self.set_atom_char(_capitalize(self.get_atom_char(idx)), idx)

        # write bonds
        for edge_list in self.graph.values():
            for edge in edge_list:
                self.set_bond_char(edge.bond_char, edge.bond_idx)


class Edge:

    def __init__(self, idx_a, idx_b, bond_char, bond_idx):
        self.idx_a = idx_a
        self.idx_b = idx_b
        self.bond_char = bond_char
        self.bond_idx = bond_idx

    def other_end(self, idx):
        if idx == self.idx_a:
            return self.idx_b
        elif idx == self.idx_b:
            return self.idx_a
        return None

    def __str__(self):  # for debugging

        return f"[{self.idx_a} -> {self.idx_b} | {self.bond_idx}]"
