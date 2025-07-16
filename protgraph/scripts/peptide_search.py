from functools import lru_cache
from math import trunc
import igraph as ig
import pandas as pd
from typing import List, Optional, Tuple

Path = Tuple[int, ...]

def contains_path_string(g: ig.Graph, peptide: str) -> list[list[int]]:
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
        #TODO: is staying in the same node needed in current architecture?
        #stay inside the same node
        while idx < L and matched < peptide_length:
            if aminoacids[idx] == peptide[matched]:
                if matched + 1 == peptide_length:   #+1 cause lists are 0-indixed
                    return (node,)  #node fully consumed, peptide found
                idx += 1
                matched += 1
            else:                                   #peptide not found
                return None

        #leave the node when its text is fully consumed
        if idx == L:
            for neighbor_node in g.neighbors(node, mode="OUT"):
                neighbor_aminoacids = g.vs[neighbor_node]["aminoacid"]
                if neighbor_aminoacids in ignore_set:
                    return None    #_end_ node has been reached, so peptide didnt match
                else:
                    if matched < peptide_length and neighbor_aminoacids[0] == peptide[matched]: #check if neighbor node is possible match
                        if matched + 1 == peptide_length: #check if already peptide found
                            return (node, neighbor_node) 
                        sub_path =  dfs(neighbor_node, 1, matched + 1) #go futher along the path
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
                    paths.append([node])  # peptide already finished   
                path = dfs(node, start + m, m)
                if path is not None:                    # peptide found 
                    paths.append(list(path))

    return paths

def get_peptides(protein_id, **kwargs):
    peptides = []
    metadata = {}
    if kwargs["peptide_file"]:
        df = pd.read_csv(kwargs["peptide_file"])
        df["Normalized Protein ID"] = df["Protein ID"].str.split("-").str[0] #sometimes the protein ids are for a specific isoform (eg. P10636-2), we only want the protein id
        normalized_isoforms = df.groupby("Normalized Protein ID")
        grouped = normalized_isoforms.get_group(protein_id)
        if kwargs["compare_columns"]:
            meta_df = pd.DataFrame()
            try:
                meta_df = pd.read_csv(kwargs["metadata_file"])
            except Exception as e:
                print("Please specify a metadata file")
            grouped = grouped.merge(meta_df, on="Sample")
            modified = pd.DataFrame()
            if kwargs["median"]:
                modified = grouped.groupby(["Sequence", kwargs["compare_columns"]])["Intensity"].median()
            else:
                modified = grouped.groupby(["Sequence", kwargs["compare_columns"]])["Intensity"].mean()
            modified = modified.unstack(level=kwargs["compare_columns"])
            peptides = modified.index
            metadata = {
                index: [value for value in row]
                for index, row in modified.iterrows()
            }
        else:
            if kwargs["median"]:
                modified = grouped.groupby(["Sequence"])["Intensity"].median()
            else:
                modified = grouped.groupby(["Sequence"])["Intensity"].mean()
            peptides = modified.index
            metadata = modified.to_dict()
    else: 
        peptides = kwargs["peptide"]
    return peptides, metadata

def metadata_to_string(metadata):
    if type(metadata) == list:
        metadata = [str(trunc(i) if trunc(i) == i else i) for i in metadata] #dont want to have .0
        metadata = ", ".join(metadata)
        metadata = "(" + metadata + ")"
    else: 
        metadata = str(trunc(metadata) if trunc(metadata) == metadata else metadata)
    return metadata

def add_peptides_to_graph(graph, peptides, metadata, show_intensity, merge_peptides, count):
    node_peptides = dict()
    edge_peptides = dict()
    for peptide in peptides:
        paths = contains_path_string(graph, peptide)     
        for path in paths:
            index = 0 #TODO: custom dict with new method for these duplicate three lines
            node = path[index]
            node_matched_peptides = node_peptides.pop(node, set())
            node_matched_peptides.add(peptide)
            node_peptides[node] = node_matched_peptides
            while index < len(path)-1:
                index += 1
                node = path[index]
                node_matched_peptides = node_peptides.pop(node, set())
                node_matched_peptides.add(peptide)
                node_peptides[node] = node_matched_peptides
                predecessor = path[index-1]
                edge_id = graph.get_eid(predecessor, node)
                edge_matched_peptides = edge_peptides.pop(edge_id, set())
                edge_matched_peptides.add(peptide)
                edge_peptides[edge_id] = edge_matched_peptides

    #add peptide information to the spanning nodes and edges
    for key, value in edge_peptides.items():
        peptides = list(value)  #converting the set to a list and sorting to have a deterministic seqeuence of the peptides 
        peptides.sort()
        if count:
            graph.es[key]["count"] = len(peptides)
        peptide_string = ", ".join(peptides)
        graph.es[key]["peptides"] = peptide_string
        if show_intensity:
            intensity = []
            for peptide in peptides:
                if peptide in metadata.keys():
                    intensity.append(metadata_to_string(metadata[peptide]))
            intensity_string = ", ".join(intensity)
            graph.es[key]["intensity"] = intensity_string
    for key, value in node_peptides.items():
        peptides = list(value)
        peptides.sort()
        if merge_peptides:
            result = []
            peptides.sort(key = len)#sort ascending by length of the peptide to allow the intelligent skip with i>j in the loop
            for i, peptide in enumerate(peptides):  #check for each peptide i if there is a longer peptide j which fully contains peptide i
                if not any(peptide in possible_super_peptide for j, possible_super_peptide in enumerate(peptides) if i < j):
                    result.append(peptide)
            peptides = result
        if count:
            graph.vs[key]["count"] = len(peptides)
        peptide_string = ", ".join(peptides)
        graph.vs[key]["peptides"] = peptide_string
        if show_intensity:
            intensity = []
            for peptide in peptides:
                if peptide in metadata.keys():
                    intensity.append(metadata_to_string(metadata[peptide]))
            intensity_string = ", ".join(intensity)
            graph.vs[key]["intensity"] = intensity_string
    return