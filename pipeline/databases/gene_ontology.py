from typing import Callable, Generator, Union

import networkx as nx

from databases import uniprot
from download import download
from enrichment import test, correction

ORGANISM = {"file": {9606: "human"}}


def get_ontology(
    namespaces: list[str] = [
        "cellular_compartment", "molecular_function", "biological_process"
    ]
) -> Generator[dict[str, Union[str, list[str]]], None, None]:
    """
    Yields Gene Ontology terms from the given namespaces.

    Args:
        namespaces: The Gene Ontology namespaces to consider terms from.

    Yields:
        A dictionary containing a Gene Ontology terms id, name, namespace and 
            related terms.
    """
    term = {}
    for line in download.txt("http://purl.obolibrary.org/obo/go.obo"):
        if any(
                line.startswith("{}:".format(tag))
                for tag in ("format-version", "data-version", "subsetdef",
                            "synonymtypedef", "default-namespace", "ontology",
                            "property_value")):
            continue
        elif line == "[Term]" or line == "[Typedef]":
            if term.get("id") and term.get("name") and term.get(
                    "namespace") in namespaces:
                yield term
            term = {}
        elif any(
                line.startswith("{}:".format(tag))
                for tag in ("id", "name", "namespace")):
            term[line.split(":",
                            maxsplit=1)[0]] = line.split(":",
                                                         maxsplit=1)[1].strip()
        elif line.startswith("is_a:"):
            if "is_a" not in term:
                term["is_a"] = []
            term["is_a"].append(
                line.split(":", maxsplit=1)[1].split("!")[0].strip())


def get_annotation(
    taxonomy_identifier: int = 9606,
    namespaces: list[str] = ["C", "F", "P"]
) -> Generator[tuple[str, str], None, None]:
    """
    Yields Gene Ontology annotations within specified namespaces.

    Args:
        taxonomy_identifier: The taxonomy identifier.
        namespace: The Gene Ontology namespace identifiers.

    Yields:
        Pairs of protein accessions and Gene Ontology term identifiers.
    """
    primary_accession = uniprot.get_primary_accession(taxonomy_identifier)

    for row in download.tabular_txt(
            "http://geneontology.org/gene-associations/goa_{organism}.gaf.gz".
            format(organism=ORGANISM["file"][taxonomy_identifier]),
            skiprows=41,
            delimiter="\t",
            usecols=[0, 1, 4, 8, 12]):
        if row[0] == "UniProtKB" and row[8] in namespaces and row[
                12] == "taxon:{}".format(taxonomy_identifier):
            for protein in primary_accession.get(row[1], {row[1]}):
                yield (protein, row[4])

    for row in download.tabular_txt(
            "http://geneontology.org/gene-associations/goa_{organism}_isoform.gaf.gz"
            .format(organism=ORGANISM["file"][taxonomy_identifier]),
            skiprows=41,
            delimiter="\t",
            usecols=[0, 4, 8, 12, 16]):
        if row[0] == "UniProtKB" and row[8] in namespaces and row[
                12] == "taxon:{}".format(
                    taxonomy_identifier) and row[16].startswith("UniProtKB:"):
            yield (row[16].split(":")[1], row[4])


def get_enrichment(
        networks: list[nx.Graph],
        test: Callable[[int, int, int, int], float] = test.hypergeometric,
        correction: Callable[[dict[int, float]],
                             dict[int, float]] = correction.benjamini_hochberg,
        taxonomy_identifier: int = 9606,
        namespaces: list[str] = [
            "cellular_compartment", "molecular_function", "biological_process"
        ],
        annotation_as_reference: bool = True
) -> dict[nx.Graph, dict[str, float]]:
    """
    Test the protein-protein interaction networks for enrichment of Gene 
    Ontology terms.

    Args:
        networks: The protein-protein interaction networks.
        test: The statistical test used to assess enrichment of a Gene Ontology 
            term.
        correction: The procedure to correct for multiple testing of multiple 
            terms and networks.
        taxonomy_identifier: The taxonomy identifier.
        namespaces: The Gene Ontology namespaces.
        annotation_as_reference: If True, compute enrichment with respect to the 
            entire species-specific Gene Ontology annotation in namespaces, else 
            with respect to the union of the protein-protein interaction 
            networks.

    Returns:
        Adjusted p-values for the enrichment of each Gene Ontology term by each 
        network.
    """
    annotation = {}
    for protein, term in get_annotation(taxonomy_identifier, [{
            "cellular_component": "C",
            "molecular_function": "F",
            "biological_process": "P"
    }[ns] for ns in namespaces]):
        if annotation_as_reference or any(
                protein in network.nodes() for network in networks):
            if term not in annotation:
                annotation[term] = set()
            annotation[term].add(protein)

    annotation = {
        term: proteins for term, proteins in annotation.items() if proteins
    }

    name = {}
    for term in get_ontology(namespaces):
        if term["id"] in annotation:
            name[term["id"]] = term["name"]

    annotated_proteins = set.union(*annotation.values())

    annotated_network_proteins = {
        network: len(annotated_proteins.intersection(network.nodes()))
        for network in networks
    }

    intersection = {
        network: {
            term: len(annotation[term].intersection(network.nodes()))
            for term in annotation
        } for network in networks
    }

    p_value = correction({(network, term):
                          test(intersection[network][term],
                               len(annotated_proteins), len(annotation[term]),
                               annotated_network_proteins[network])
                          for term in annotation for network in networks})

    return {
        network: {(term, name[term]): p_value[(network, term)]
                  for term in annotation} for network in networks
    }
