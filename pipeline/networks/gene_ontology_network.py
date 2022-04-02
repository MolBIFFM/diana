"""
Gene Ontology network

Nodes are Gene Ontology terms annotated with proteins from a protein-protein
interaction network species. Edges are directed term relationships within the
Gene Ontology.
"""

from typing import Callable, Container

import networkx as nx
from analysis import correction, test
from databases import gene_ontology


def get_gene_ontology_network(
        protein_protein_interaction_network: nx.Graph,
        namespaces: Container[str] = ("cellular_component",
                                      "molecular_function",
                                      "biological_process"),
        enrichment_test: Callable[[int, int, int, int],
                                  float] = test.hypergeometric,
        multiple_testing_correction: Callable[[dict[str, float]], dict[
            str, float]] = correction.benjamini_hochberg,
        taxonomy_identifier: int = 9606) -> nx.Graph:
    """
    Assemble a Gene Ontology network corresponding to the protein-protein
    interaction network.

    Args:
        protein_protein_interaction_network: The protein-protein interaction
            network.
        namespaces: The Gene Ontology namespaces.
        enrichment_test: The statistical test used to assess enrichment of a
            term by the protein-protein interaction network.
        multiple_testing_correction: The procedure to correct for testing of
            multiple terms.
        taxonomy_identifier: The taxonomy identifier.

    Returns:
        The Gene Ontology network.
    """
    network, go_id = nx.DiGraph(), {}
    for term in gene_ontology.get_ontology(namespaces):
        network.add_node(term["id"])
        network.nodes[term["id"]]["term"] = term["name"]
        network.nodes[term["id"]]["namespace"] = term["namespace"]

        for parent in term.get("is_a", []):
            network.add_edge(term["id"], parent)

        for alt_id in term["alt_id"]:
            if alt_id not in go_id:
                go_id[alt_id] = set()
            go_id[alt_id].add(term["id"])

    annotation = {}
    for protein, term in gene_ontology.get_annotation(
            taxonomy_identifier, gene_ontology.convert_namespaces(namespaces)):
        for primary_term in go_id.get(term, {term}):
            if primary_term not in annotation:
                annotation[primary_term] = set()
            annotation[primary_term].add(protein)

    network.remove_nodes_from(
        [term for term in network if term not in annotation])

    annotated_proteins = set.union(*annotation.values())

    network_intersection = {
        term: len(annotation[term].intersection(
            protein_protein_interaction_network.nodes())) for term in annotation
    }

    p_value = multiple_testing_correction({
        term: enrichment_test(
            network_intersection[term], len(annotated_proteins),
            len(annotation[term]),
            len(
                annotated_proteins.intersection(
                    protein_protein_interaction_network.nodes())))
        for term in network
    })

    network.remove_nodes_from(
        [term for term in annotation if not network_intersection[term]])

    for term in network:
        network.nodes[term]["proteins"] = network_intersection[term]
        network.nodes[term]["p-value"] = p_value[term]

    return network


def get_term_sizes(network: nx.Graph) -> dict[str, int]:
    """
    Returns the sizes of Gene Ontology term annotation.

    Args:
        network: The Gene Ontology network.

    Returns:
        The number of proteins associated with any term in the Gene Ontology
        network.
    """
    return {term: network.nodes[term]["proteins"] for term in network}


def export(network: nx.Graph, basename: str) -> None:
    """
    Exports the Gene Ontology network.

    Args:
        network: The Gene Ontology network.
        basename: The base file name.
    """
    nx.write_graphml_xml(network,
                         f"{basename}.graphml",
                         named_key_ids=True,
                         infer_numeric_types=True)