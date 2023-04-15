"""
Gene Ontology network

Nodes are Gene Ontology terms associated with proteins from a species of
interest. Edges are directed term relationships within the Gene Ontology.
"""
import os
from typing import (Callable, Collection, Container, Hashable, Iterable,
                    Literal, Mapping, Optional)

import networkx as nx
import scipy.stats
from algorithms import correction
from databases import gene_ontology


def get_network(proteins: Iterable[str],
                reference: Optional[Container[str]] = None,
                namespaces: Collection[
                    Literal["cellular_component", "molecular_function",
                            "biological_process"]] = frozenset(),
                enrichment_test: Callable[[int, int, int, int], float] = lambda
                k, M, n, N: float(scipy.stats.hypergeom.sf(k - 1, M, n, N)),
                multiple_testing_correction: Callable[
                    [dict[Hashable, float]],
                    Mapping[Hashable, float]] = correction.benjamini_yekutieli,
                organism: int = 9606,
                file_ontology: Optional[str] = None,
                file_annotation: Optional[str] = None,
                file_annotation_isoform: Optional[str] = None,
                file_uniprot: Optional[str] = None) -> nx.DiGraph:
    """
    Assemble a Gene Ontology network from proteins.

    Args:
        proteins: The proteins to assemble the Gene Ontology network from.
        reference: Optional reference set of proteins with respect to which
            enrichment is computed. If not provided, the entire Gene Ontology
            annotation specific to the organism of interest is used.
        namespaces: The Gene Ontology namespaces to consider. If empty, any
            namespace is considered.
        enrichment_test: The statistical test used to assess enrichment of a
            term by the protein-protein interaction network.
        multiple_testing_correction: The procedure to correct for testing of
            multiple terms.
        organism: The NCBI taxonomy identifier for the organism of interest.
        file_ontology: The optional local file location to parse terms from.
        file_annotation: The optional local file location to parse annotations
            from.
        file_annotation: The optional local file location to parse isoform
            annotations from.
        file_uniprot: The optional local file location to parse accessions from.

    Returns:
        The Gene Ontology network.
    """
    # Initialize the Gene Ontology network.
    network = nx.DiGraph()

    # Add Gene Ontology terms and their relations to the Gene Ontology network
    # and compile a map from alternative to primary Gene Ontology term
    # identifiers.
    go_id: dict[str, set[str]] = {}
    for term in gene_ontology.get_ontology(namespaces, file_ontology):
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

    # Compile a map from primary Gene Ontology term identifiers to UniProt
    # accessions of annotated proteins.
    annotations: dict[str, set[str]] = {}
    for protein, annotated_term in gene_ontology.get_annotation(
            organism, gene_ontology.convert_namespaces(namespaces),
            file_annotation, file_annotation_isoform, file_uniprot):
        for primary_term in go_id.get(annotated_term, {annotated_term}):
            if primary_term not in annotations:
                annotations[primary_term] = set()

            if reference is None or protein in reference:
                annotations[primary_term].add(protein)

    # Remove nodes representing Gene Ontology terms not associated with any
    # proteins.
    network.remove_nodes_from([
        term for term in network
        if term not in annotations or not annotations[term]
    ])

    # Discard Gene Ontology terms not associated with any proteins.
    annotations = {
        term: annotation
        for term, annotation in annotations.items()
        if annotation
    }

    # Compile a reference set of proteins annotated with any Gene Ontology term.
    annotated_proteins = set.union(*annotations.values())

    # Compile a map from individual Gene Ontology terms to subsets of associated
    # proteins.
    prt_intersection = {
        term: annotation.intersection(proteins)
        for term, annotation in annotations.items()
    }

    # Compute p-values for the enrichment of individual Gene Ontology terms
    # corrected for testing multiple Gene Ontology terms.
    p_value = multiple_testing_correction({
        term: enrichment_test(len(prt_intersection[term]),
                              len(annotated_proteins), len(annotations[term]),
                              len(annotated_proteins.intersection(proteins)))
        for term in network
    })

    # Annotate the nodes of the Gene Ontology network with associated proteins.
    for node in network:
        network.nodes[node]["p-value"] = p_value[node]
        network.nodes[node]["number of associated proteins"] = len(
            prt_intersection[node])
        network.nodes[node]["associated proteins"] = " ".join(
            sorted(prt_intersection[node]))

    # Return the Gene Ontology network.
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
    # Return a map from Gene Ontology terms to the number of associated
    # proteins.
    return {
        term: network.nodes[term]["number of associated proteins"]
        for term in network
    }


def export(network: nx.Graph, basename: str) -> Optional[str]:
    """
    Exports the Gene Ontology network if possible without overwriting an
    existing file.


    Args:
        network: The Gene Ontology network.
        basename: The base file name.

    Returns:
        The file the Gene Ontology network was exported to if there is no naming
        conflict.
    """
    # Avoid overwriting an existing file.
    if os.path.isfile(f"{basename}.graphml"):
        return None

    # Export the element tree of the Gene Ontology network.
    nx.write_graphml_xml(network,
                         f"{basename}.graphml",
                         named_key_ids=True,
                         infer_numeric_types=True)

    # Return the file name of the Gene Ontology network.
    return f"{basename}.graphml"
