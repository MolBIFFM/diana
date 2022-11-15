"""
Gene Ontology network

Nodes are Gene Ontology terms associated with proteins from a species of
interest. Edges are directed term relationships within the Gene Ontology.
"""

from typing import (Callable, Collection, Container, Hashable, Iterable,
                    Mapping, Optional)

import networkx as nx
import scipy.stats
from algorithms import correction
from databases import gene_ontology


def get_network(proteins: Iterable[str],
                reference: Optional[Container[str]] = None,
                namespaces: Collection[str] = ("cellular_component",
                                               "molecular_function",
                                               "biological_process"),
                enrichment_test: Callable[[int, int, int, int], float] = lambda
                k, M, n, N: scipy.stats.hypergeom.sf(k - 1, M, n, N),
                multiple_testing_correction: Callable[
                    [dict[Hashable, float]],
                    Mapping[Hashable, float]] = correction.benjamini_yekutieli,
                organism: int = 9606) -> nx.DiGraph:
    """
    Assemble a Gene Ontology network from proteins.

    Args:
        proteins: The proteins to assemble the Gene Ontology network from.
        reference: Optional reference set of proteins with respect to which
            enrichment is computed. If not provided, the entire Gene Ontology
            annotation specific to the organism of interest is used.
        namespaces: The Gene Ontology namespaces.
        enrichment_test: The statistical test used to assess enrichment of a
            term by the protein-protein interaction network.
        multiple_testing_correction: The procedure to correct for testing of
            multiple terms.
        organism: The NCBI taxonomy identifier for the organism of interest.

    Returns:
        The Gene Ontology network.
    """
    network = nx.DiGraph()
    go_id: dict[str, set[str]] = {}
    for term in gene_ontology.get_ontology(namespaces):
        if isinstance(term["id"], str) and isinstance(term["namespace"], str):
            network.add_node(term["id"])
            network.nodes[term["id"]]["term"] = term["name"]
            network.nodes[term["id"]]["namespace"] = term["namespace"].replace(
                "_", " ")

            for parent in term.get("is_a", []):
                network.add_edge(term["id"], parent)

            for alt_id in term["alt_id"]:
                if alt_id not in go_id:
                    go_id[alt_id] = set()
                go_id[alt_id].add(term["id"])

    annotation: dict[str, set[str]] = {}
    for protein, annotated_term in gene_ontology.get_annotation(
            organism, gene_ontology.convert_namespaces(namespaces)):
        for primary_term in go_id.get(annotated_term, {annotated_term}):
            if primary_term not in annotation:
                annotation[primary_term] = set()

            if reference is None or protein in reference:
                annotation[primary_term].add(protein)

    network.remove_nodes_from(
        [term for term in network if term not in annotation])

    annotated_proteins = set.union(*annotation.values())

    network_intersection = {
        term: term_proteins.intersection(proteins)
        for term, term_proteins in annotation.items()
    }

    p_value = multiple_testing_correction({
        term: enrichment_test(len(network_intersection[term]),
                              len(annotated_proteins), len(annotation[term]),
                              len(annotated_proteins.intersection(proteins)))
        for term in network
    })

    for node in network:
        network.nodes[node]["p-value"] = p_value[node]
        network.nodes[node]["number of proteins"] = len(
            network_intersection[node])
        network.nodes[node]["proteins"] = " ".join(
            sorted(network_intersection[node]))

    return network


def get_term_sizes(network: nx.Graph) -> dict[str, int]:
    """
    Returns the sizes of Gene Ontology term annotation.

    Args:
        network: The Gene Ontology network.

    Returns:
        The number of proteins from the initial protein-protein interaction
        network associated with any term in the Gene Ontology
        network.
    """
    return {term: network.nodes[term]["number of proteins"] for term in network}


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
