"""
Reactome network

Nodes are Reactome pathways associated with proteins from a species of interest.
Edges are directed pathway relationships within Reactome.
"""
import os
from typing import Callable, Container, Hashable, Iterable, Mapping, Optional

import networkx as nx
import scipy.stats
from algorithms import correction
from databases import reactome


def get_network(proteins: Iterable[str],
                reference: Optional[Container[str]] = None,
                enrichment_test: Callable[[int, int, int, int], float] = lambda
                k, M, n, N: float(scipy.stats.hypergeom.sf(k - 1, M, n, N)),
                multiple_testing_correction: Callable[
                    [dict[Hashable, float]],
                    Mapping[Hashable, float]] = correction.benjamini_yekutieli,
                organism: int = 9606,
                file_pathways: Optional[str] = None,
                file_pathways_relation: Optional[str] = None,
                file_accession_map: Optional[str] = None,
                file_uniprot: Optional[str] = None) -> nx.DiGraph:
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
        file_pathways: The optional local file location to parse pathways
            from.
        file_pathways_relation: The optional local file location to parse
            pathway relations from.
        file_accession_map: The optional local file location to parse accession
            associations from.
        file_uniprot: The optional local file location to parse accessions from.

    Returns:
        The Reactome network.
    """
    # Initialize the Reactome network.
    network = nx.DiGraph()

    # Add Reactome pathways to the Reactome network.
    for pathway, name in reactome.get_pathways(organism, file_pathways):
        network.add_node(pathway)
        network.nodes[pathway]["pathway"] = name

    # Add Reactome pathway relations to the Reactome network.
    for parent, child in reactome.get_pathway_relations(file_pathways_relation):
        if parent in network and child in network:
            network.add_edge(child, parent)

    # Compile a map from Reactome pathway identifiers to UniProt accessions of
    # annotated proteins.
    annotations: dict[str, set[str]] = {}
    for protein, pathway in reactome.get_pathway_annotation(
            organism, file_accession_map, file_uniprot):
        if pathway not in annotations:
            annotations[pathway] = set()

        if not reference or protein in reference:
            annotations[pathway].add(protein)

    # Remove nodes representing Reactome pathways not associated with any
    # proteins.
    network.remove_nodes_from([
        pathway for pathway in network
        if pathway not in annotations or not annotations[pathway]
    ])

    # Discard Reactome pathways not associated with any proteins.
    annotations = {
        pathway: annotation
        for pathway, annotation in annotations.items()
        if annotation
    }

    # Compile a reference set of proteins annotated with any Reactome pathway.
    annotated_proteins = set.union(*annotations.values())

    # Compile a map from individual Reactome pathways to subsets of associated
    # proteins.
    prt_intersection = {
        pathway: annotation.intersection(proteins)
        for pathway, annotation in annotations.items()
    }

    # Compute p-values for the enrichment of individual Reactome pathways
    # corrected for testing multiple Reactome pathways.
    p_value = multiple_testing_correction({
        pathway: enrichment_test(len(prt_intersection[pathway]),
                                 len(annotated_proteins),
                                 len(annotations[pathway]),
                                 len(annotated_proteins.intersection(proteins)))
        for pathway in network
    })

    # Annotate the nodes of the Reactome network with associated proteins.
    for node in network:
        network.nodes[node]["p-value"] = p_value[node]
        network.nodes[node]["number of associated proteins"] = len(
            prt_intersection[node])
        network.nodes[node]["associated proteins"] = " ".join(
            sorted(prt_intersection[node]))

    # Return the Reactome network.
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
    # Return a map from Reactome pathways to the number of associated proteins.
    return {
        pathway: network.nodes[pathway]["number of associated proteins"]
        for pathway in network
    }


def export(network: nx.Graph, basename: str) -> Optional[str]:
    """
    Exports the Reactome network if possible without overwriting an existing
    file.


    Args:
        network: The Reactome network.
        basename: The base file name.

    Returns:
        The file the Reactome network was exported to if there is no naming
        conflict.
    """
    # Avoid overwriting an existing file.
    if os.path.isfile(f"{basename}.graphml"):
        return None

    # Export the element tree of the Reactome network.
    nx.write_graphml_xml(network,
                         f"{basename}.graphml",
                         named_key_ids=True,
                         infer_numeric_types=True)

    # Return the file name of the Reactome network.
    return f"{basename}.graphml"
