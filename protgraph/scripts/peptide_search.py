from functools import lru_cache
import igraph as ig
from typing import List, Optional, Tuple

Path = Tuple[int, ...]

def contains_path_string(g: ig.Graph, peptide: str) -> List[List[int]]:
    """Returns a list all paths where a DAG path concatenates to `peptide` under the start/end cut rule."""
    peptide_length = len(peptide)
    ignore_set = {"__start__", "__end__"}
    paths = []

    @lru_cache(maxsize=None)
    def dfs(node: int, idx: int, matched: int) -> Optional[Path]:
        """
        Return True if from state (node node, cursor idx inside its text,
        matched chars of search) we can finish the pattern.
        """
        aminoacids = g.vs[node]["aminoacid"]

        #_end_ node has been reached, so peptide didnt match
        if aminoacids in ignore_set:
            return None

        L = len(aminoacids)

        #1) stay inside the same node
        while idx < L and matched < peptide_length:
            if aminoacids[idx] == peptide[matched]:
                if matched + 1 == peptide_length:   #+1 cause lists are 0-indixed
                    #TODO: add check if node fully consumed
                    return (node,)                     #peptide found
                idx += 1
                matched += 1
            else:                                   #peptide not found
                return None

        #2) leave the node when its text is fully consumed
        if idx == L:
            for neighbor_node in g.neighbors(node, mode="OUT"):
                neighbor_aminoacids = g.vs[neighbor_node]["aminoacid"]
                if neighbor_aminoacids in ignore_set:
                    return None    #_end_ node has been reached, so peptide didnt match
                else:
                    if matched < peptide_length and neighbor_aminoacids[0] == peptide[matched]: #check if neighbor node is possible match
                        if matched + 1 == peptide_length:
                            return (node, neighbor_node) #peptide found
                        sub_path =  dfs(neighbor_node, 1, matched + 1)
                        if sub_path is not None:
                            return (node,) + sub_path
        return None

    #seed all possible start positions
    for node in range(g.vcount()):
        aminoacids = g.vs[node]["aminoacid"]
        if aminoacids in ignore_set:
            continue

        L = len(aminoacids)
        # any suffix of aminoacids that matches a prefix of peptide may be a start
        for start in range(L):
            m = 0
            while m < peptide_length and start + m < L and aminoacids[start + m] == peptide[m]:
                m += 1
            if m:
                if m == peptide_length:
                    paths.append([node])                # peptide already finished  TODO: from when to where     
                path = dfs(node, start + m, m)
                if path is not None:                    # peptide found TODO: from when
                    paths.append([node] + list(path))

    return paths