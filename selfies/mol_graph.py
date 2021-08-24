import functools
from typing import Dict, List, Optional, Union

from selfies.bond_constraints import get_bonding_capacity
from selfies.kekulize import find_perfect_matching, prune_to_pi_subgraph


class Atom:

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


class MolecularDFSTree:

    def __init__(self):
        self._atoms = list()
        self._bond_dict = dict()
        self._adj_list = list()
        self._bond_counts = list()
        self._ring_bond_flags = list()
        self._aro_subgraph = dict()

    def __len__(self):
        return len(self._atoms)

    def has_bond(self, a: int, b: int) -> bool:
        return ((a, b) in self._bond_dict) or ((b, a) in self._bond_dict)

    def has_out_ring_bond(self, src: int) -> bool:
        return self._ring_bond_flags[src]

    def get_atom(self, idx: int) -> Atom:
        return self._atoms[idx]

    def get_atoms(self) -> List[Atom]:
        return self._atoms

    def get_dirbond(self, src, dst) -> DirectedBond:
        return self._bond_dict[(src, dst)]

    def get_bond_count(self, idx: int) -> int:
        return self._bond_counts[idx]

    def get_out_dirbonds(self, src: int) -> List[DirectedBond]:
        return self._adj_list[src]

    def get_aromatic_subgraph(self) -> Dict[int, List[int]]:
        return self._aro_subgraph

    def add_atom(self, atom: Atom) -> None:
        atom.index = len(self)
        self._atoms.append(atom)
        self._adj_list.append(list())
        self._bond_counts.append(0)
        self._ring_bond_flags.append(False)
        if atom.is_aromatic:
            self._aro_subgraph[atom.index] = list()

    def add_bond(
            self, src: int, dst: int,
            order: Union[int, float], stereo: str
    ) -> None:
        assert src < dst

        bond = DirectedBond(src, dst, order, stereo, False)
        self._add_bond_at_loc(bond, -1)
        self._bond_counts[src] += order
        self._bond_counts[dst] += order

        if order == 1.5:
            self._aro_subgraph.setdefault(src, []).append(dst)
            self._aro_subgraph.setdefault(dst, []).append(src)

    def add_placeholder(self, src: int) -> int:
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
            self._aro_subgraph.setdefault(a, []).append(b)
            self._aro_subgraph.setdefault(b, []).append(a)

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

    def is_kekulized(self) -> bool:
        return not bool(self._aro_subgraph)

    def kekulize(self) -> bool:
        if self.is_kekulized():
            return True

        # prune to pi subgraph
        pi_subgraph = prune_to_pi_subgraph(self)
        self._aro_subgraph = dict()

        # cleaning up bond counts
        for i in range(len(self)):
            self._bond_counts[i] = int(self._bond_counts[i])

        # find matching and make double bonds
        matching = find_perfect_matching(pi_subgraph)
        if matching is None:
            return False

        for (a, b) in matching:
            self.update_bond_order(a, b, 2)
        return True

    def _add_bond_at_loc(self, bond, pos):
        self._bond_dict[(bond.src, bond.dst)] = bond

        out_edges = self._adj_list[bond.src]
        if (pos == -1) or (pos == len(out_edges)):
            out_edges.append(bond)
        elif out_edges[pos] is None:
            out_edges[pos] = bond
        else:
            out_edges.insert(pos, bond)
