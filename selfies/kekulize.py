from typing import Dict, Iterable, List, Set, Tuple, Union

from selfies.grammar_rules import find_element, get_num_from_bond, \
    parse_atom_symbol

ATOM_TYPE = 1
BRANCH_TYPE = 2
RING_TYPE = 3


def kekulize_parser(smiles_gen: Iterable[Tuple[str, str, int]]) \
        -> Iterable[Tuple[str, str, int]]:
    """Kekulizes a SMILES in the form of an iterable.

    This method intercepts the output of ``encoder._parse_smiles``, and
    acts as filter that kekulizes the SMILES. The motivation for having
    this setup is that string parsing and concatenation is minimized,
    as the parsing is already done by ``_parse_smiles``.

    Reference: https://depth-first.com/articles/2020/02/10/a-comprehensive
               -treatment-of-aromaticity-in-the-smiles-language/

    :param smiles_gen: an iterator returned by ``encoder._parse_smiles``.
    :return: an iterator representing the kekulized SMILES, in the same
        format as that returned by ``encoder._parse_smiles``.
    """

    # save to list, so the iterator can be used across multiple functions
    # change elements from tuple -> list to allow in-place modifications
    smiles_symbols = list(map(list, smiles_gen))

    mol_graph = MolecularGraph(smiles_symbols)

    rings = {}
    _build_molecular_graph(mol_graph, smiles_symbols, rings)

    if mol_graph.aro_indices:
        _kekulize(mol_graph)

    for x in mol_graph.smiles_symbols:  # return as iterator
        yield tuple(x)


def _build_molecular_graph(graph,
                           smiles_symbols: List[List[Union[str, int]]],
                           rings: Dict[int, Tuple[int, int]],
                           prev_idx: int = -1,
                           curr_idx: int = -1) -> int:
    """From the iterator returned by ``encoder._parse_smiles``, builds
    a graph representation of the molecule.

    This is done by iterating through ``smiles_symbols``, and then adding bonds
    to the molecular graph. Note that ``smiles_symbols`` is mutated in this
    method, for convenience.

    :param graph: the MolecularGraph to be added to.
    :param smiles_symbols: a list created from the iterator returned
        by ``encoder._parse_smiles``.
    :param rings: an, initially, empty dictionary used to keep track of
        rings to be made.
    :param prev_idx:
    :param curr_idx:
    :return: the last index of ``smiles_symbols`` that was processed.
    """

    while curr_idx + 1 < len(smiles_symbols):

        curr_idx += 1
        _, symbol, symbol_type = smiles_symbols[curr_idx]

        if symbol_type == ATOM_TYPE:
            if prev_idx >= 0:
                graph.add_bond(prev_idx, curr_idx, curr_idx)
            prev_idx = curr_idx

        elif symbol_type == BRANCH_TYPE:
            if symbol == '(':
                curr_idx = _build_molecular_graph(graph, smiles_symbols, rings,
                                                  prev_idx, curr_idx)
            else:
                break

        else:
            if symbol in rings:
                left_idx, left_bond_idx = rings.pop(symbol)
                right_idx, right_bond_idx = prev_idx, curr_idx

                # we mutate one bond index to be '', so that we
                # can faithfully represent the bond to be localized at
                # one index. For example, C=1CCCC=1 --> C1CCCC=1.

                if smiles_symbols[left_bond_idx][0] != '':
                    bond_idx = left_bond_idx
                    smiles_symbols[right_bond_idx][0] = ''
                else:
                    bond_idx = right_bond_idx
                    smiles_symbols[left_bond_idx][0] = ''

                graph.add_bond(left_idx, right_idx, bond_idx)
            else:
                rings[symbol] = (prev_idx, curr_idx)

    return curr_idx


def _kekulize(mol_graph) -> None:
    """Kekulizes the molecular graph.

    :param mol_graph: a molecular graph to be kekulized.
    :return: None.
    """

    mol_graph.prune_to_pi_subgraph()

    visited = set()
    for i in mol_graph.get_nodes_by_num_edges():
        success = mol_graph.dfs_assign_bonds(i, visited, set(), set())
        if not success:
            raise ValueError("Kekulization Failed.")

    mol_graph.write_to_smiles_symbols()


# Aromatic Helper Methods and Classes

# key = aromatic SMILES element, value = number of valence electrons
# Note: wild card '*' not supported currently
_aromatic_valences = {
    'b': 3, 'c': 4, 'n': 5, 'p': 5, 'as': 5,
    'o': 6, 's': 6, 'se': 6
}


def _capitalize(atom_symbol: str) -> str:
    """Capitalizes the element portion of an aromatic SMILES atom symbol,
    converting it into a standard SMILES atom symbol.

    :param atom_symbol: an aromatic SMILES atom symbol.
    :return: the capitalized ``atom_symbol``.
    """

    s, _ = find_element(atom_symbol)
    return atom_symbol[:s] + atom_symbol[s].upper() + atom_symbol[s + 1:]


def _is_aromatic(atom_symbol: str) -> bool:
    """Checks whether a SMILES atom symbol is an aromatic SMILES atom symbol.

    An aromatic SMILES atom symbol is indicated by an element substring
    that is not capitalized.

    :param atom_symbol: a SMILES atom symbol.
    :return: True, if ``atom_symbol`` is an aromatic atom symbol,
        and False otherwise.
    """

    s, e = find_element(atom_symbol)

    if e == len(atom_symbol):  # optimization to prevent string copying
        element = atom_symbol
    else:
        element = atom_symbol[s: e]

    if element[0].isupper():  # check if element is capitalized
        return False

    if element not in _aromatic_valences:
        raise ValueError("Kekulization Failed: aromatic symbol {} "
                         "not recognized.".format(atom_symbol))
    return True


def _in_pi_subgraph(atom_symbol: str, bonds: Tuple[str]) -> bool:
    """Checks whether a SMILES atom symbol should be a node in the pi
    subgraph, based on its bonds.

    More specifically, an atom should be a node in the pi subgraph if it has
    an unpaired valence electron, and thus, is able to make a double bond.

    Reference: https://depth-first.com/articles/2020/02/10/a-comprehensive
               -treatment-of-aromaticity-in-the-smiles-language/

    :param atom_symbol: a SMILES atom symbol representing an atom.
    :param bonds: the bonds connected to ``atom_symbol``.
    :return: True if ``atom_symbol`` should be included in the pi subgraph,
        and False otherwise.
    """

    atom, h_count, charge = parse_atom_symbol(atom_symbol)

    used_electrons = 0
    for b in bonds:
        used_electrons += get_num_from_bond(b)

    # e.g. c1ccccc1
    if (atom == 'c') and (h_count == charge == 0) and (len(bonds) == 2):
        h_count += 1  # implied bonded hydrogen

    if h_count > 1:
        raise ValueError("Kekulization Failed: "
                         "{} not supported.".format(atom_symbol))

    elif h_count == 1:  # e.g. [nH]
        used_electrons += 1

    valence = _aromatic_valences[atom] - charge
    free_electrons = valence - used_electrons
    return free_electrons % 2 != 0


class MolecularGraph:
    """A molecular graph.

    This molecular graph operates based on the ``smiles_symbols`` data
    structure. Indices from this list represent nodes or edges, depending
    on whether they point to a SMILES atom(s) or bond.

    :ivar smiles_symbols: the list created from the iterator returned by
        ``encoder._parse_smiles``. Serves as the base data structure
        of this class, as everything is communicated through indices
        referring to elements of this list.
    :ivar graph: the key is an index of the atom(s) from ``smiles_symbols``.
        The value is a list of Bond objects representing the connected
        bonds. Represents the actual molecular graph.
    :ivar aro_indices: a set of indices of atom(s) from ``smiles_symbols``
        that are aromatic in the molecular graph.
    """

    def __init__(self, smiles_symbols: List[List[Union[str, int]]]):
        self.smiles_symbols = smiles_symbols
        self.graph = {}
        self.aro_indices = set()

    def get_atom_symbol(self, idx: int) -> str:
        """Getter that returns the SMILES symbol representing an atom
        at a specified index.

        :param idx: an index in ``smiles_symbols``.
        :return: the SMILES symbol representing an atom at index
            ``idx`` in ``smiles_symbols``.
        """

        return self.smiles_symbols[idx][1]

    def get_bond_symbol(self, idx: int) -> str:
        """Getter that returns the SMILES symbol representing a bond at
        a specified index.

        :param idx: an index in ``smiles_symbols``.
        :return: the SMILES symbol representing a bond at index
            ``idx`` in ``smiles_symbols``.
        """

        return self.smiles_symbols[idx][0]

    def get_nodes_by_num_edges(self) -> List[int]:
        """Returns all nodes (or indices) stored in this molecular graph
        in a semi-sorted order by number of edges.

        This is to optimize the speed of ``dfs_assign_bonds``; starting
        with nodes that have fewer edges will improve computational time
        as there are fewer bond configurations to explore. Instead of fully
        sorting the returned list, a compromise is made, and nodes with exactly
        one edge are added to the list's beginning.

        :return: a list of the nodes (or indices) of this molecular graph,
            semi-sorted by number of edges.
        """

        ends = []  # nodes with exactly 1 edge
        middles = []  # nodes with 2+ edges

        for idx, edges in self.graph.items():
            if len(edges) > 1:
                middles.append(idx)
            else:
                ends.append(idx)

        ends.extend(middles)
        return ends

    def set_atom_symbol(self, atom_symbol: str, idx: int) -> None:
        """Setter that updates the SMILES symbol representing an atom(s) at
        a specified index.

        :param atom_symbol: the new value of the atom symbol at ``idx``.
        :param idx: an index in ``smiles_symbols``.
        :return: None.
        """

        self.smiles_symbols[idx][1] = atom_symbol

    def set_bond_symbol(self, bond_symbol: str, idx: int) -> None:
        """Setter that updates the SMILES symbol representing a bond at
        a specified index.

        :param bond_symbol: the new value of the bond symbol at ``idx``.
        :param idx: an index in ``smiles_symbols``.
        :return: None.
        """

        self.smiles_symbols[idx][0] = bond_symbol

    def add_bond(self, idx_a: int, idx_b: int, bond_idx: int) -> None:
        """Adds a bond (or edge) to this molecular graph between atoms
        (or nodes) at two specified indices.

        :param idx_a: the index of one atom (or node) of this bond.
        :param idx_b:the index of one atom (or node) of this bond.
        :param bond_idx: the index of this bond.
        :return: None.
        """

        atom_a = self.get_atom_symbol(idx_a)
        atom_b = self.get_atom_symbol(idx_b)
        atom_a_aro = (idx_a in self.aro_indices) or _is_aromatic(atom_a)
        atom_b_aro = (idx_b in self.aro_indices) or _is_aromatic(atom_b)
        bond_symbol = self.get_bond_symbol(bond_idx)

        if atom_a_aro:
            self.aro_indices.add(idx_a)

        if atom_b_aro:
            self.aro_indices.add(idx_b)

        if bond_symbol == ':':
            self.aro_indices.add(idx_a)
            self.aro_indices.add(idx_b)

            # Note: ':' bonds are edited here to ''
            self.set_bond_symbol('', bond_idx)
            bond_symbol = ''

        edge = Bond(idx_a, idx_b, bond_symbol, bond_idx)

        self.graph.setdefault(idx_a, []).append(edge)
        self.graph.setdefault(idx_b, []).append(edge)

    def prune_to_pi_subgraph(self) -> None:
        """Removes nodes and edges from this molecular graph such that
        it becomes the pi subgraph.

        The remaining graph will only contain aromatic atoms (or nodes)
        that belong in the pi-subgraph, and the bonds that are aromatic
        and between such atoms.

        :return: None.
        """

        # remove non-aromatic nodes
        non_aromatic = self.graph.keys() - self.aro_indices
        for i in non_aromatic:
            self.graph.pop(i)

        # remove non-pi subgraph nodes
        for i in self.aro_indices:

            atom = self.get_atom_symbol(i)
            bonds = tuple(edge.bond_symbol for edge in self.graph[i])

            if not _in_pi_subgraph(atom, bonds):
                self.graph.pop(i)

        # remove irrelevant edges
        for idx, edges in self.graph.items():

            keep = list(filter(lambda e: (e.idx_a in self.graph)
                                         and (e.idx_b in self.graph)
                                         and (e.bond_symbol == ''),
                               edges))
            self.graph[idx] = keep

    def dfs_assign_bonds(self, idx: int,
                         visited: Set[int],
                         matched_nodes: Set[int],
                         matched_edges) -> bool:
        """After calling ``prune_to_pi_subgraph``, this method assigns
        double bonds between pairs of nodes such that every node is
        paired or matched.

        This is done recursively in a depth-first search fashion.

        :param idx: the index of the current atom (or node).
        :param visited: a set of the indices of nodes that have been visited.
        :param matched_nodes: a set of the indices of nodes that have been
            matched, i.e., assigned a double bond.
        :param matched_edges: a set of the bonds that have been matched.
        :return: True, if a valid bond assignment was found; False otherwise.
        """

        if idx in visited:
            return True

        edges = self.graph[idx]

        if idx in matched_nodes:

            # recursively try to match adjacent nodes. If the matching
            # fails, then we must backtrack.
            visited_save = visited.copy()

            visited.add(idx)
            for e in edges:
                adj = e.other_end(idx)
                if not self.dfs_assign_bonds(adj, visited,
                                             matched_nodes,
                                             matched_edges):
                    visited &= visited_save
                    return False
            return True

        else:

            # list of candidate edges that can become a double bond
            candidates = list(
                filter(lambda i: i.other_end(idx) not in matched_nodes, edges)
            )

            if not candidates:
                return False  # idx is unmatched, but all adj nodes are matched

            matched_edges_save = matched_edges.copy()

            for e in candidates:

                # match nodes connected by c
                matched_nodes.add(e.idx_a)
                matched_nodes.add(e.idx_b)
                matched_edges.add(e)

                success = self.dfs_assign_bonds(idx, visited,
                                                matched_nodes,
                                                matched_edges)

                if success:
                    e.bond_symbol = '='
                    return True
                else:  # the matching failed, so we must backtrack

                    for edge in matched_edges - matched_edges_save:
                        edge.bond_symbol = ''
                        matched_nodes.discard(edge.idx_a)
                        matched_nodes.discard(edge.idx_b)

                    matched_edges &= matched_edges_save

            return False

    def write_to_smiles_symbols(self):
        """Updates and mutates ``self.smiles_symbols`` with the information
         contained in ``self.graph``.

        After kekulizing the molecular graph, this method is called to
        merge the new information back into the original data structure.

        :return: None.
        """

        # capitalize aromatic molecules
        for idx in self.aro_indices:
            self.set_atom_symbol(_capitalize(self.get_atom_symbol(idx)), idx)

        # write bonds
        for edge_list in self.graph.values():
            for edge in edge_list:
                bond_symbol = edge.bond_symbol
                bond_idx = edge.bond_idx

                self.set_bond_symbol(bond_symbol, bond_idx)

                # branches record the next symbol as their bond, so we
                # must update accordingly
                if (bond_idx > 0) and \
                        (self.smiles_symbols[bond_idx - 1][2] == BRANCH_TYPE):
                    self.set_bond_symbol(bond_symbol, bond_idx - 1)


class Bond:
    """Represents a bond or edge in MolecularGraph.

    Recall that the following indices are with respect to ``smiles_symbols``
    in MolecularGraph.

    :ivar idx_a: the index of one atom or node of this bond.
    :ivar idx_b: the index of one atom or node of this bond.
    :ivar bond_symbol: the SMILES symbol representing this bond (e.g. '#').
    :ivar bond_idx: the index of this bond or edge.
    """

    def __init__(self, idx_a, idx_b, bond_symbol, bond_idx):
        self.idx_a = idx_a
        self.idx_b = idx_b
        self.bond_symbol = bond_symbol
        self.bond_idx = bond_idx

    def __eq__(self, other):
        if isinstance(other, type(self)):
            return (self.idx_a, self.idx_b) == (other.idx_a, other.idx_b)
        return NotImplemented

    def __hash__(self):
        return hash((self.idx_a, self.idx_b))

    def other_end(self, idx):
        """Given an index representing one end of this bond, returns
        the index representing the other end.

        :param idx: an index of one atom or node of this bond.
        :return: the index of the other atom or node of this bond, or
            None if ``idx`` is an invalid input.
        """

        if idx == self.idx_a:
            return self.idx_b
        elif idx == self.idx_b:
            return self.idx_a
        return None
