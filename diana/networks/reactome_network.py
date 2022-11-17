"""
Reactome network

Nodes are Reactome pathways associated with proteins from a species of interest.
Edges are directed pathway relationships within Reactome.
"""

from typing import Callable, Container, Hashable, Iterable, Mapping, Optional

import networkx as nx
import scipy.stats
from algorithms import correction
from databases import reactome


def get_network(proteins: Iterable[str],
                reference: Optional[Container[str]] = None,
                enrichment_test: Callable[[int, int, int, int], float] = lambda
                k, M, n, N: scipy.stats.hypergeom.sf(k - 1, M, n, N),
                multiple_testing_correction: Callable[
                    [dict[Hashable, float]],
                    Mapping[Hashable, float]] = correction.benjamini_yekutieli,
                organism: int = 9606) -> nx.DiGraph:
    """
    Assemble a Reactome network from proteins.

    Args:
        proteins: The proteins to assemble the Reactome network from.
        reference: Optional reference set of proteins with respect to which
            enrichment is computed. If not provided, the entire Reactome pathway
            annotation specific to the organism of interest is used.
        enrichment_test: The statistical test used to assess enrichment of a
            pathway by the protein-protein interaction network.
        multiple_testing_correction: The procedure to correct for testing of
            multiple pathways.
        organism: The NCBI taxonomy identifier for the organism of interest.

    Returns:
        The Reactome network.
    """
    network = nx.DiGraph()
    for pathway, name in reactome.get_pathways(organism):
        network.add_node(pathway)
        network.nodes[pathway]["pathway"] = name

    for child, parent in reactome.get_pathway_relations():
        if child in network and parent in network:
            network.add_edge(child, parent)

    annotation: dict[str, set[str]] = {}
    for protein, pathway in reactome.get_pathway_annotation(organism):
        if pathway not in annotation:
            annotation[pathway] = set()

        if not reference or protein in reference:
            annotation[pathway].add(protein)

    network.remove_nodes_from([
        pathway for pathway in network
        if pathway not in annotation or not annotation[pathway]
    ])

    annotation = {
        pathway: annot for pathway, annot in annotation.items() if annot
    }

    annotated_proteins = set.union(*annotation.values())

    network_intersection = {
        pathway: annot.intersection(proteins)
        for pathway, annot in annotation.items()
    }

    p_value = multiple_testing_correction({
        pathway: enrichment_test(len(network_intersection[pathway]),
                                 len(annotated_proteins),
                                 len(annotation[pathway]),
                                 len(annotated_proteins.intersection(proteins)))
        for pathway in network
    })

    for node in network:
        network.nodes[node]["p-value"] = p_value[node]
        network.nodes[node]["number of proteins"] = len(
            network_intersection[node])
        network.nodes[node]["proteins"] = " ".join(
            sorted(network_intersection[node]))

    return network


def get_pathway_sizes(network: nx.Graph) -> dict[str, int]:
    """
    Returns the sizes of Reactome pathway annotation.

    Args:
        network: The Reactome network.

    Returns:
        The number of proteins from the initial protein-protein interaction
        network associated with any pathway in Reactome.
    """
    return {
        pathway: network.nodes[pathway]["number of proteins"]
        for pathway in network
    }


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
