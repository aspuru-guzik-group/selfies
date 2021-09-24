import heapq
import itertools
from collections import deque
from typing import List, Optional


def find_perfect_matching(graph: List[List[int]]) -> Optional[List[int]]:
    """Finds a perfect matching for an undirected graph (without self-loops).

    :param graph: an adjacency list representing the input graph.
    :return: a list representing a perfect matching, where j is the i-th
        element if nodes i and j are matched. Returns None, if the graph cannot
        be perfectly matched.
    """

    # start with a maximal matching for efficiency
    matching = _greedy_matching(graph)

    unmatched = set(i for i in range(len(graph)) if matching[i] is None)
    while unmatched:

        # find augmenting path which starts at root
        root = unmatched.pop()
        path = _find_augmenting_path(graph, root, matching)

        if path is None:
            return None
        else:
            _flip_augmenting_path(matching, path)
            unmatched.discard(path[0])
            unmatched.discard(path[-1])

    return matching


def _greedy_matching(graph):
    matching = [None] * len(graph)
    free_degrees = [len(graph[i]) for i in range(len(graph))]
    # free_degrees[i] = number of unmatched neighbors for node i

    # prioritize nodes with fewer unmatched neighbors
    node_pqueue = [(free_degrees[i], i) for i in range(len(graph))]
    heapq.heapify(node_pqueue)

    while node_pqueue:
        _, node = heapq.heappop(node_pqueue)

        if (matching[node] is not None) or (free_degrees[node] == 0):
            continue  # node cannot be matched

        # match node with first unmatched neighbor
        mate = next(i for i in graph[node] if matching[i] is None)
        matching[node] = mate
        matching[mate] = node

        for adj in itertools.chain(graph[node], graph[mate]):
            free_degrees[adj] -= 1
            if (matching[adj] is None) and (free_degrees[adj] > 0):
                heapq.heappush(node_pqueue, (free_degrees[adj], adj))

    return matching


def _find_augmenting_path(graph, root, matching):
    assert matching[root] is None

    # run modified BFS to find path from root to unmatched node
    other_end = None
    node_queue = deque([root])

    # parent BFS tree - None indicates an unvisited node
    parents = [None] * len(graph)
    parents[root] = [None, None]

    while node_queue:
        node = node_queue.popleft()

        for adj in graph[node]:
            if matching[adj] is None:  # unmatched node
                if adj != root:  # augmenting path found!
                    parents[adj] = [node, adj]
                    other_end = adj
                    break
            else:
                adj_mate = matching[adj]
                if parents[adj_mate] is None:  # adj_mate not visited
                    parents[adj_mate] = [node, adj]
                    node_queue.append(adj_mate)

        if other_end is not None:
            break  # augmenting path found!

    if other_end is None:
        return None
    else:
        path = []
        node = other_end
        while node != root:
            path.append(parents[node][1])
            path.append(parents[node][0])
            node = parents[node][0]
        return path


def _flip_augmenting_path(matching, path):
    for i in range(0, len(path), 2):
        a, b = path[i], path[i + 1]
        matching[a] = b
        matching[b] = a
