"""
Reactome network

Nodes are Reactome pathways annotated with proteins from a protein-protein
interaction network. Edges are directed pathway relationships within Reactome.
"""

from typing import Callable

import networkx as nx
from analysis import correction, test
from databases import reactome


def get_reactome_network(protein_protein_interaction_network: nx.Graph,
                         enrichment_test: Callable[[int, int, int, int],
                                                   float] = test.hypergeometric,
                         multiple_testing_correction: Callable[
                             [dict[str, float]],
                             dict[str, float]] = correction.benjamini_hochberg,
                         taxonomy_identifier: int = 9606) -> nx.Graph:
    """
    Assemble a Reactome network corresponding to the protein-protein interaction
    network.

    Args:
        protein_protein_interaction_network: The protein-protein interaction
            network.
        enrichment_test: The statistical test used to assess enrichment of a
            pathway by the protein-protein interaction network.
        multiple_testing_correction: The procedure to correct for testing of
            multiple pathways.
        taxonomy_identifier: The taxonomy identifier.

    Returns:
        The Reactome network.
    """
    network = nx.DiGraph()
    for pathway, name in reactome.get_pathways(taxonomy_identifier):
        network.add_node(pathway)
        network.nodes[pathway]["pathway"] = name

    for child, parent in reactome.get_pathway_relations():
        if child in network and parent in network:
            network.add_edge(child, parent)

    pathways = {}
    for protein, pathway in reactome.get_pathway_annotation(
            taxonomy_identifier):
        if pathway not in pathways:
            pathways[pathway] = set()
        pathways[pathway].add(protein)

    network.remove_nodes_from(
        [pathway for pathway in network if pathway not in pathways])

    annotated_proteins = set.union(*pathways.values())

    annotated_network_proteins = {
        pathway: len(pathways[pathway].intersection(
            protein_protein_interaction_network.nodes()))
        for pathway in pathways
    }

    p_value = multiple_testing_correction({
        pathway: enrichment_test(
            annotated_network_proteins[pathway], len(annotated_proteins),
            len(pathways[pathway]),
            len(
                annotated_proteins.intersection(
                    protein_protein_interaction_network.nodes())))
        for pathway in network
    })

    network.remove_nodes_from([
        pathway for pathway in pathways
        if not annotated_network_proteins[pathway]
    ])

    for pathway in network:
        network.nodes[pathway]["proteins"] = annotated_network_proteins[pathway]
        network.nodes[pathway]["p-value"] = p_value[pathway]

    return network


def get_pathway_sizes(network: nx.Graph) -> dict[str, int]:
    """
    Returns the sizes of Reactome pathway annotation.

    Args:
        network: The Reactome network.

    Returns:
        The number of proteins associated with any pathway in Reactome.
    """
    return {pathway: network.nodes[pathway]["proteins"] for pathway in network}


def export(network: nx.Graph, basename: str) -> None:
    """
    Exports the Reactome network.

    Args:
        network: The Reactome network.
        basename: The base file name.
    """
    nx.write_graphml_xml(network,
                         f"{basename}.graphml",
                         named_key_ids=True,
                         infer_numeric_types=True)
