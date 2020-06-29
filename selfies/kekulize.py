from __future__ import annotations

from typing import Dict, Iterable, List, Set, Tuple, Union

from selfies.grammar_rules import get_num_from_bond

ATOM_TYPE = 1
BRANCH_TYPE = 2
RING_TYPE = 3


def kekulize_parser(smiles_gen: Iterable[str, str, int]) \
        -> Iterable[str, str, int]:
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
    smiles_chars = list(map(list, smiles_gen))

    mol_graph = MolecularGraph(smiles_chars)

    rings = {}
    _build_molecular_graph(mol_graph, smiles_chars, rings)

    if mol_graph.aro_indices:
        _kekulize(mol_graph)

    for x in mol_graph.smiles_chars:  # return as iterator
        yield tuple(x)


def _build_molecular_graph(graph: MolecularGraph,
                           smiles_chars: List[List[Union[str, int]]],
                           rings: Dict[int, Tuple[int, int]],
                           prev_idx: int = -1,
                           curr_idx: int = -1) -> int:
    """From the iterator returned by ``encoder._parse_smiles``, builds
    a graph representation of the molecule.

    This is done by iterating through ``smiles_chars``, and then adding bonds
    to the molecular graph. Note that ``smiles_chars`` is mutated in this
    method, for convenience.

    :param graph: the MolecularGraph to be added to.
    :param smiles_chars: a list created from the iterator returned
        by ``encoder._parse_smiles``.
    :param rings: an, initially, empty dictionary used to keep track of
        rings to be made.
    :param prev_idx:
    :param curr_idx:
    :return: the last index of ``smiles_char`` that was processed.
    """

    while curr_idx + 1 < len(smiles_chars):

        curr_idx += 1
        _, char, char_type = smiles_chars[curr_idx]

        if char_type == ATOM_TYPE:
            if prev_idx >= 0:
                graph.add_bond(prev_idx, curr_idx, curr_idx)
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

                # we mutate one bond index to be '', so that we
                # can faithfully represent the bond to be localized at
                # one index. For example, C=1CCCC=1 --> C1CCCC=1.

                if smiles_chars[left_bond_idx][0] != '':
                    bond_idx = left_bond_idx
                    smiles_chars[right_bond_idx][0] = ''
                else:
                    bond_idx = right_bond_idx
                    smiles_chars[left_bond_idx][0] = ''

                graph.add_bond(left_idx, right_idx, bond_idx)
            else:
                rings[char] = (prev_idx, curr_idx)

    return curr_idx


def _kekulize(mol_graph: MolecularGraph) -> None:
    """Kekulizes the molecular graph.

    :param mol_graph: a molecular graph to be kekulized.
    :return: None.
    """

    mol_graph.prune_to_pi_subgraph()

    visited = set()
    for i in mol_graph.get_nodes_by_num_edges():
        matched = set()

        success = mol_graph.dfs_assign_bonds(i, visited, matched)
        if not success:
            raise ValueError("Kekulization Failed.")

    mol_graph.write_to_smiles_chars()


# Aromatic Helper Methods and Classes

# key = aromatic SMILES element, value = number of valence electrons
# Note: wild card '*' not supported currently
_aromatic_valences = {
    'b': 3, 'c': 4, 'n': 5, 'p': 5, 'as': 5,
    'o': 6, 's': 6, 'se': 6
}


def _find_element(char: str) -> Tuple[int, int]:
    """Returns the indices of the element component of a SMILES character.

    That is, if char[i:j] is the element substring of the SMILES character,
    then (i, j) is returned. For example:
        *   _find_element('b') = (0, 1).
        *   _find_element('B') = (0, 1).
        *   _find_element('[13C]') = (3, 4).
        *   _find_element('[nH+]') = (1, 2).

    :param char: a SMILES character.
    :return: a tuple of the indices of the element substring of ``char``.
    """

    if char[0] != '[':
        return 0, len(char)

    i = 1
    while char[i].isdigit():  # skip isotope number
        i += 1

    if char[i + 1].isalpha() and char[i + 1] != 'H':
        return i, i + 2
    else:
        return i, i + 1


def _parse_char(char: str) -> Tuple[str, int, int]:
    """Parses a SMILES character and returns its element component,
    hydrogen count, and charge.

    See http://opensmiles.org/opensmiles.html for the formal grammar
    of SMILES characters. Note that only @ and @@ are currently supported
    as chiral symbols.

    :param char: a SMILES character.
    :return: a tuple of (1) the element of ``char``, (2) the hydrogen count,
        and (3) the charge.
    """

    if char[0] != '[':
        return char, 0, 0

    atom_start, atom_end = _find_element(char)
    i = atom_end

    # skip chirality
    if char[i] == '@':  # e.g. @
        i += 1
    if char[i] == '@':  # e.g. @@
        i += 1

    h_count = 0  # hydrogen count
    if char[i] == 'H':
        h_count = 1

        i += 1
        if char[i].isdigit():  # e.g. [CH2]
            h_count = int(char[i])
            i += 1

    charge = 0  # charge count
    if char[i] in ('+', '-'):
        charge = 1 if char[i] == '+' else -1

        i += 1
        if char[i] in ('+', '-'):  # e.g. [Cu++]
            while char[i] in ('+', '-'):
                charge += (1 if char[i] == '+' else -1)
                i += 1

        elif char[i].isdigit():  # e.g. [Cu+2]
            s = i
            while char[i].isdigit():
                i += 1
            charge *= int(char[s:i])

    return char[atom_start: atom_end], h_count, charge


def _capitalize(char: str) -> str:
    """Capitalizes the element portion of an aromatic SMILES character,
    converting it into a standard SMILES character.

    :param char: an aromatic SMILES character.
    :return: the capitalized ``char``.
    """

    c, _ = _find_element(char)
    return char[:c] + char[c].upper() + char[c + 1:]


def _is_aromatic(char: str) -> bool:
    """Checks whether a SMILES character is an aromatic SMILES character.

    An aromatic SMILES character is indicated by an element substring
    that is not capitalized.

    :param char: a SMILES character.
    :return: True, if ``char`` is an aromatic character, and False otherwise.
    """

    s, e = _find_element(char)

    if e == len(char):  # optimization to prevent string copying
        element = char
    else:
        element = char[s: e]

    if element[0].isupper():  # check if element is capitalized
        return False

    if element not in _aromatic_valences:
        raise ValueError(f"Kekulization Failed: aromatic symbol {char} "
                         f"not recognized.")
    return True


def _in_pi_subgraph(char: str, bonds: Tuple[str]) -> bool:
    """Checks whether a SMILES character should be a node in the pi
    subgraph, based on its bonds.

    More specifically, an atom should be a node in the pi subgraph if it has
    an unpaired valence electron, and thus, is able to make a double bond.

    Reference: https://depth-first.com/articles/2020/02/10/a-comprehensive
               -treatment-of-aromaticity-in-the-smiles-language/

    :param char: a SMILES character representing an atom.
    :param bonds: the bonds connected to ``char``.
    :return: True if ``char`` should be included in the pi subgraph,
        and False otherwise.
    """

    atom, h_count, charge = _parse_char(char)

    used_electrons = 0
    for b in bonds:
        used_electrons += get_num_from_bond(b)

    if char == 'c' and len(bonds) == 2:  # e.g. c1ccccc1
        h_count += 1  # implied bonded hydrogen

    if h_count > 1:
        raise ValueError(f"Kekulization Failed: {char} not supported.")

    elif h_count == 1:  # e.g. [nH]
        used_electrons += 1

    valence = _aromatic_valences[atom] - charge
    free_electrons = valence - used_electrons
    return free_electrons % 2 != 0


class MolecularGraph:
    """A molecular graph.

    This molecular graph operates based on the ``smiles_chars`` data
    structure. Indices from this list represent nodes or edges, depending
    on whether they point to a SMILES atom(s) or bond.

    :ivar smiles_chars: the list created from the iterator returned by
        ``encoder._parse_smiles``. Serves as the base data structure
        of this class, as everything is communicated through indices
        referring to elements of this list.
    :ivar graph: the key is an index of the atom(s) from ``smiles_chars``.
        The value is a list of Bond objects representing the connected
        bonds. Represents the actual molecular graph.
    :ivar aro_indices: a set of indices of atom(s) from ``smiles_chars``
        that are aromatic in the molecular graph.
    """
    smiles_chars: List[List[Union[str, int]]]
    graph: Dict[int, List[Bond]]
    aro_indices: Set[int]

    def __init__(self, smiles_chars: List[List[Union[str, int]]]):
        self.smiles_chars = smiles_chars
        self.graph = {}
        self.aro_indices = set()

    def get_atom_char(self, idx: int) -> str:
        """Getter that returns the SMILES character representing an atom(s)
        at a specified index.

        :param idx: an index in ``smiles_chars``.
        :return: the SMILES character representing an atom(s) at index
            ``idx`` in ``smiles_chars``.
        """

        return self.smiles_chars[idx][1]

    def get_bond_char(self, idx: int) -> str:
        """Getter that returns the SMILES character representing a bond at
        a specified index.

        :param idx: an index in ``smiles_chars``.
        :return: the SMILES character representing a bond at index
            ``idx`` in ``smiles_chars``.
        """

        return self.smiles_chars[idx][0]

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

    def set_atom_char(self, atom_char: str, idx: int) -> None:
        """Setter that updates the SMILES character representing an atom(s) at
        a specified index.

        :param atom_char: the new value of the atom character at ``idx``.
        :param idx: an index in ``smiles_chars``.
        :return: None.
        """

        self.smiles_chars[idx][1] = atom_char

    def set_bond_char(self, bond_char: str, idx: int) -> None:
        """Setter that updates the SMILES character representing a bond at
        a specified index.

        :param bond_char: the new value of the bond character at ``idx``.
        :param idx: an index in ``smiles_chars``.
        :return: None.
        """

        self.smiles_chars[idx][0] = bond_char

    def add_bond(self, idx_a: int, idx_b: int, bond_idx: int) -> None:
        """Adds a bond (or edge) to this molecular graph between atoms
        (or nodes) at two specified indices.

        :param idx_a: the index of one atom (or node) of this bond.
        :param idx_b:the index of one atom (or node) of this bond.
        :param bond_idx: the index of this bond.
        :return: None.
        """

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

            # Note: ':' bonds are edited here to ''
            self.set_bond_char('', bond_idx)
            bond_char = ''

        edge = Bond(idx_a, idx_b, bond_char, bond_idx)

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

    def dfs_assign_bonds(self, idx: int,
                         visited: Set[int],
                         matched: Set[int]) -> bool:
        """After calling ``prune_to_pi_subgraph``, this method assigns
        double bonds between pairs of nodes such that every node is
        paired or matched.

        This is done recursively in a depth-first search fashion.

        :param idx: the index of the current atom (or node).
        :param visited: a set of the indices of nodes that have been visited.
        :param matched: a set of the indices of nodes that have been matched,
            i.e., assigned a double bond.
        :return: True, if a valid bond assignment was found; False otherwise.
        """

        if idx in visited:
            return True

        edges = self.graph[idx]

        if idx in matched:

            # recursively try to match adjacent nodes. If the matching
            # fails, then we must backtrack.

            visited.add(idx)

            for e in edges:
                adj = e.other_end(idx)
                if not self.dfs_assign_bonds(adj, visited, matched):
                    visited.remove(idx)
                    return False
            return True

        else:

            # list of candidate edges that can become a double bond
            candidates = list(filter(lambda i: i.other_end(idx) not in matched,
                                     edges))

            if not candidates:
                return False  # idx is unmatched, but all adj nodes are matched

            for c in candidates:

                # match nodes connected by c
                c.bond_char = '='
                matched.add(c.idx_a)
                matched.add(c.idx_b)

                success = self.dfs_assign_bonds(idx, visited, matched)

                if success:
                    return True
                else:  # the matching failed, so we must backtrack
                    c.bond_char = ''
                    matched.remove(c.idx_a)
                    matched.remove(c.idx_b)

            return False

    def write_to_smiles_chars(self):
        """Updates and mutates ``self.smiles_chars`` with the information
         contained in ``self.graph``.

        After kekulizing the molecular graph, this method is called to
        merge the new information back into the original data structure.

        :return: None.
        """

        # capitalize aromatic molecules
        for idx in self.aro_indices:
            self.set_atom_char(_capitalize(self.get_atom_char(idx)), idx)

        # write bonds
        for edge_list in self.graph.values():
            for edge in edge_list:
                bond_char = edge.bond_char
                bond_idx = edge.bond_idx

                self.set_bond_char(bond_char, bond_idx)

                # branches record the next character as their bond, so we
                # must update accordingly
                if (bond_idx > 0) and \
                        (self.smiles_chars[bond_idx - 1][2] == BRANCH_TYPE):
                    self.set_bond_char(bond_char, bond_idx - 1)


class Bond:
    """Represents a bond or edge in MolecularGraph.

    Recall that the following indices are with respect to ``smiles_chars``
    in MolecularGraph.

    :ivar idx_a: the index of one atom or node of this bond.
    :ivar idx_b: the index of one atom or node of this bond.
    :ivar bond_char: the SMILES character representing this bond (e.g. '#').
    :ivar bond_idx: the index of this bond or edge.
    """
    idx_a: int
    idx_b: int
    bond_char: str
    bond_idx: int

    def __init__(self, idx_a, idx_b, bond_char, bond_idx):
        self.idx_a = idx_a
        self.idx_b = idx_b
        self.bond_char = bond_char
        self.bond_idx = bond_idx

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
