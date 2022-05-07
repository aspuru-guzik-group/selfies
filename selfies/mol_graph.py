import functools
import itertools
from typing import List, Optional, Union
from dataclasses import dataclass, field

from selfies.bond_constraints import get_bonding_capacity
from selfies.constants import AROMATIC_VALENCES
from selfies.utils.matching_utils import find_perfect_matching


@dataclass
class Attribution:
    """A dataclass that contains token string and its index.
    """
    #: token index
    index: int
    #: token string
    token: str


@dataclass
class AttributionMap:
    """A mapping from input to single output token showing which
    input tokens created the output token.
    """
    #: Index of output token
    index: int
    #: Output token
    token: str
    #: List of input tokens that created the output token
    attribution: List[Attribution] = field(default_factory=list)


class Atom:
    """An atom with associated specifications (e.g. charge, chirality).
    """

    def __init__(
            self,
            element: str,
            is_aromatic: bool,
            isotope: Optional[int] = None,
            chirality: Optional[str] = None,
            h_count: Optional[int] = None,
            charge: int = 0
    ):
        self.index = None
        self.element = element
        self.is_aromatic = is_aromatic
        self.isotope = isotope
        self.chirality = chirality
        self.h_count = h_count
        self.charge = charge

    @property
    @functools.lru_cache()
    def bonding_capacity(self):
        bond_cap = get_bonding_capacity(self.element, self.charge)
        bond_cap -= 0 if (self.h_count is None) else self.h_count
        return bond_cap

    def invert_chirality(self) -> None:
        if self.chirality == "@":
            self.chirality = "@@"
        elif self.chirality == "@@":
            self.chirality = "@"


class DirectedBond:
    """A bond that contains directional information.
    """

    def __init__(
            self,
            src: int,
            dst: int,
            order: Union[int, float],
            stereo: Optional[str],
            ring_bond: bool
    ):
        self.src = src
        self.dst = dst
        self.order = order
        self.stereo = stereo
        self.ring_bond = ring_bond


class MolecularGraph:
    """A molecular graph.

    Molecules can be viewed as weighted undirected graphs. However, SMILES
    and SELFIES strings are more naturally represented as weighted directed
    graphs, where the direction of the edges specifies the order of atoms
    and bonds in the string.
    """

    def __init__(self, attributable=False):
        self._roots = list()  # stores root atoms, where traversal begins
        self._atoms = list()  # stores atoms in this graph
        self._bond_dict = dict()  # stores all bonds in this graph
        self._adj_list = list()  # adjacency list, representing this graph
        self._bond_counts = list()  # stores number of bonds an atom has made
        self._ring_bond_flags = list()  # stores if an atom makes a ring bond
        self._delocal_subgraph = dict()  # delocalization subgraph
        self._attribution = dict()  # attribution of each atom/bond
        self._attributable = attributable

    def __len__(self):
        return len(self._atoms)

    def has_bond(self, a: int, b: int) -> bool:
        if a > b:
            a, b = b, a
        return (a, b) in self._bond_dict

    def has_out_ring_bond(self, src: int) -> bool:
        return self._ring_bond_flags[src]

    def get_attribution(
        self,
        o: Union[DirectedBond, Atom]
    ) -> List[Attribution]:
        if self._attributable and o in self._attribution:
            return self._attribution[o]
        return None

    def get_roots(self) -> List[int]:
        return self._roots

    def get_atom(self, idx: int) -> Atom:
        return self._atoms[idx]

    def get_atoms(self) -> List[Atom]:
        return self._atoms

    def get_dirbond(self, src, dst) -> DirectedBond:
        return self._bond_dict[(src, dst)]

    def get_out_dirbonds(self, src: int) -> List[DirectedBond]:
        return self._adj_list[src]

    def get_bond_count(self, idx: int) -> int:
        return self._bond_counts[idx]

    def add_atom(self, atom: Atom, mark_root: bool = False) -> Atom:
        atom.index = len(self)

        if mark_root:
            self._roots.append(atom.index)
        self._atoms.append(atom)
        self._adj_list.append(list())
        self._bond_counts.append(0)
        self._ring_bond_flags.append(False)
        if atom.is_aromatic:
            self._delocal_subgraph[atom.index] = list()
        return atom

    def add_attribution(
            self,
            o: Union[DirectedBond, Atom],
            attr: List[Attribution]
    ) -> None:
        if self._attributable:
            if o in self._attribution:
                self._attribution[o].extend(attr)
            else:
                self._attribution[o] = attr

    def add_bond(
            self, src: int, dst: int,
            order: Union[int, float], stereo: str
    ) -> DirectedBond:
        assert src < dst

        bond = DirectedBond(src, dst, order, stereo, False)
        self._add_bond_at_loc(bond, -1)
        self._bond_counts[src] += order
        self._bond_counts[dst] += order

        if order == 1.5:
            self._delocal_subgraph.setdefault(src, []).append(dst)
            self._delocal_subgraph.setdefault(dst, []).append(src)
        return bond

    def add_placeholder_bond(self, src: int) -> int:
        out_edges = self._adj_list[src]
        out_edges.append(None)
        return len(out_edges) - 1

    def add_ring_bond(
            self, a: int, b: int,
            order: Union[int, float],
            a_stereo: Optional[str], b_stereo: Optional[str],
            a_pos: int = -1, b_pos: int = -1
    ) -> None:
        a_bond = DirectedBond(a, b, order, a_stereo, True)
        b_bond = DirectedBond(b, a, order, b_stereo, True)
        self._add_bond_at_loc(a_bond, a_pos)
        self._add_bond_at_loc(b_bond, b_pos)
        self._bond_counts[a] += order
        self._bond_counts[b] += order
        self._ring_bond_flags[a] = True
        self._ring_bond_flags[b] = True

        if order == 1.5:
            self._delocal_subgraph.setdefault(a, []).append(b)
            self._delocal_subgraph.setdefault(b, []).append(a)

    def update_bond_order(
            self, a: int, b: int,
            new_order: Union[int, float]
    ) -> None:
        assert 1 <= new_order <= 3

        if a > b:
            a, b = b, a  # swap so that a < b
        a_to_b = self._bond_dict[(a, b)]  # prev step guarantees existence
        if new_order == a_to_b.order:
            return
        elif a_to_b.ring_bond:
            b_to_a = self._bond_dict[(b, a)]
            bonds = (a_to_b, b_to_a)
        else:
            bonds = (a_to_b,)

        old_order = bonds[0].order
        for bond in bonds:
            bond.order = new_order
        self._bond_counts[a] += (new_order - old_order)
        self._bond_counts[b] += (new_order - old_order)

    def _add_bond_at_loc(self, bond, pos):
        self._bond_dict[(bond.src, bond.dst)] = bond

        out_edges = self._adj_list[bond.src]
        if (pos == -1) or (pos == len(out_edges)):
            out_edges.append(bond)
        elif out_edges[pos] is None:
            out_edges[pos] = bond
        else:
            out_edges.insert(pos, bond)

    def is_kekulized(self) -> bool:
        return not self._delocal_subgraph

    def kekulize(self) -> bool:
        # Algorithm based on Depth-First article by Richard L. Apodaca
        # Reference:
        #   https://depth-first.com/articles/2020/02/10/
        #   a-comprehensive-treatment-of-aromaticity-in-the-smiles-language/

        if self.is_kekulized():
            return True

        ds = self._delocal_subgraph
        kept_nodes = set(itertools.filterfalse(self._prune_from_ds, ds))

        # relabel kept DS nodes to be 0, 1, 2, ...
        label_to_node = list(sorted(kept_nodes))
        node_to_label = {v: i for i, v in enumerate(label_to_node)}

        # pruned and relabelled DS
        pruned_ds = [list() for _ in range(len(kept_nodes))]
        for node in kept_nodes:
            label = node_to_label[node]
            for adj in filter(lambda v: v in kept_nodes, ds[node]):
                pruned_ds[label].append(node_to_label[adj])

        matching = find_perfect_matching(pruned_ds)
        if matching is None:
            return False

        # de-aromatize and then make double bonds
        for node in ds:
            for adj in ds[node]:
                self.update_bond_order(node, adj, new_order=1)
            self._atoms[node].is_aromatic = False
            self._bond_counts[node] = int(self._bond_counts[node])

        for matched_labels in enumerate(matching):
            matched_nodes = tuple(label_to_node[i] for i in matched_labels)
            self.update_bond_order(*matched_nodes, new_order=2)

        self._delocal_subgraph = dict()  # clear DS
        return True

    def _prune_from_ds(self, node):
        adj_nodes = self._delocal_subgraph[node]
        if not adj_nodes:
            return True  # aromatic atom with no aromatic bonds

        atom = self._atoms[node]
        valences = AROMATIC_VALENCES[atom.element]

        # each bond in DS has order 1.5 - we treat them as single bonds
        used_electrons = int(self._bond_counts[node] - 0.5 * len(adj_nodes))

        if atom.h_count is None:  # account for implicit Hs
            assert atom.charge == 0
            return any(used_electrons == v for v in valences)
        else:
            valence = valences[-1] - atom.charge
            used_electrons += atom.h_count
            free_electrons = valence - used_electrons
            return not ((free_electrons >= 0) and (free_electrons % 2 != 0))
