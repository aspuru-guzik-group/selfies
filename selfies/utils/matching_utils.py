from typing import Optional, Set, Tuple


def find_perfect_matching(graph) -> Optional[Set[Tuple[int, int]]]:
    # sort nodes in graph
    bins = [[], []]
    for idx, edges in graph.items():
        if len(edges) > 1:
            bins[1].append(idx)
        else:
            bins[0].append(idx)
    sorted_nodes = bins[0] + bins[1]

    # run dfs
    visited = set()
    matching = set()
    for idx in sorted_nodes:
        if not _dfs_match_nodes(graph, idx, visited, set(), matching):
            return None

    return matching


def _dfs_match_nodes(graph, curr, visited, matched, matching):
    if curr in visited:
        return True

    neighbors = graph[curr]
    if curr in matched:
        visited_save = visited.copy()

        visited.add(curr)
        for adj in neighbors:
            if not _dfs_match_nodes(graph, adj, visited, matched, matching):
                visited &= visited_save
                return False
        return True

    else:
        candidates = list(filter(lambda x: x not in matched, neighbors))

        if not candidates:
            return False

        matching_save = matching.copy()
        for adj in candidates:
            edge = (min(curr, adj), max(curr, adj))
            matched.add(curr)
            matched.add(adj)
            matching.add(edge)

            if _dfs_match_nodes(graph, curr, visited, matched, matching):
                return True
            else:
                for edge in matching - matching_save:
                    matched.discard(edge[0])
                    matched.discard(edge[1])
                    matching.discard(edge)

        return False