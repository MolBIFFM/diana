"""The interface for the Gene Ontology."""
from typing import Generator, Union

import networkx as nx

from databases import uniprot
from download import download
from analysis import test, correction

ORGANISM = {"file": {9606: "human"}}


def get_ontology(
    namespaces: list[str] = [
        "cellular_component", "molecular_function", "biological_process"
    ]
) -> Generator[dict[str, Union[str, list[str]]], None, None]:
    """
    Yields Gene Ontology terms from the given namespaces.

    Args:
        namespaces: The Gene Ontology namespaces to consider terms from.

    Yields:
        Mappings containing a Gene Ontology terms' GO ID, name, namespace, 
            related terms and alternative GO IDs.
    """
    term = {"id": "", "name": "", "namespace": "", "is_a": [], "alt_id": []}
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
            term = {
                "id": "",
                "name": "",
                "namespace": "",
                "is_a": [],
                "alt_id": []
            }
        elif any(
                line.startswith("{}:".format(tag))
                for tag in ("id", "name", "namespace")):
            term[line.split(":",
                            maxsplit=1)[0]] = line.split(":",
                                                         maxsplit=1)[1].strip()
        elif line.startswith("is_a:"):
            term["is_a"].append(
                line.split(":", maxsplit=1)[1].split("!")[0].strip())

        elif line.startswith("alt_id:"):
            term["alt_id"].append(line.split(":", maxsplit=1)[1].strip())


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
        if row[0] == "UniProtKB" and row[3] in namespaces and row[4].split(
                ":")[-1] == str(taxonomy_identifier):
            for protein in primary_accession.get(row[1], {row[1]}):
                yield (protein, row[2])

    for row in download.tabular_txt(
            "http://geneontology.org/gene-associations/goa_{organism}_isoform.gaf.gz"
            .format(organism=ORGANISM["file"][taxonomy_identifier]),
            skiprows=41,
            delimiter="\t",
            usecols=[0, 4, 8, 12, 16]):
        if row[0] == "UniProtKB" and row[2] in namespaces and row[3].split(
                ":")[-1] == str(taxonomy_identifier) and row[4].startswith(
                    "UniProtKB:"):
            yield (row[4].split(":")[1], row[1])


def convert_namespaces(
    namespaces: list[str] = [
        "cellular_component", "molecular_function", "biological_process"
    ]
) -> list[str]:
    """
    Converts Gene Ontology namespace identifiers.

    Args:
        namespaces: Gene Ontology namespaces.

    Returns:
        The corresponding identifiers used in annotation files.
    """
    return [{
        "cellular_component": "C",
        "molecular_function": "F",
        "biological_process": "P"
    }[ns] for ns in namespaces]
