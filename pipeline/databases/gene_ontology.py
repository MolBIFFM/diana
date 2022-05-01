"""The interface for the Gene Ontology database."""
from typing import Container, Iterator, Union

from access import iterate

from databases import uniprot

ORGANISM = {"files": {9606: "human"}}


def get_ontology(namespaces: Container[str] = (
    "cellular_component", "molecular_function",
    "biological_process")) -> Iterator[dict[str, Union[str, tuple[str]]]]:
    """
    Yields Gene Ontology terms from the given namespaces.

    Args:
        namespaces: The Gene Ontology namespaces to consider terms from.

    Yields:
        Mappings containing a Gene Ontology terms' GO ID, name, namespace,
            related terms and alternative GO IDs.
    """
    term = {"id": "", "name": "", "namespace": "", "is_a": [], "alt_id": []}
    for line in iterate.txt("http://purl.obolibrary.org/obo/go.obo"):
        if any(
                line.startswith(f"{tag}:")
                for tag in ("format-version", "data-version", "subsetdef",
                            "synonymtypedef", "default-namespace", "ontology",
                            "property_value")):
            continue
        elif line in ("[Term]", "[Typedef]"):
            if term.get("id") and term.get("namespace") in namespaces:
                for attribute in ("is_a", "alt_id"):
                    term[attribute] = tuple(term[attribute])

                yield term

            term = {
                "id": "",
                "name": "",
                "namespace": "",
                "is_a": [],
                "alt_id": []
            }

        elif any(
                line.startswith(f"{tag}:")
                for tag in ("id", "name", "namespace")):
            term[line.split(":")[0]] = line.split(":", maxsplit=1)[1].strip()
        elif line.startswith("is_a:"):
            term["is_a"].append(
                line.split(":", maxsplit=1)[1].split("!")[0].strip())

        elif line.startswith("alt_id:"):
            term["alt_id"].append(line.split(":", maxsplit=1)[1].strip())


def get_annotation(
    organism: int = 9606,
    namespaces: Container[str] = ("C", "F", "P")
) -> Iterator[tuple[str, str]]:
    """
    Yields Gene Ontology annotations within specified namespaces.

    Args:
        organism: The NCBI taxonomy identifier for the organism of interest. 
        namespace: The Gene Ontology namespace identifiers.

    Yields:
        Pairs of protein accessions and Gene Ontology term identifiers.
    """
    primary_accession = uniprot.get_primary_accession(organism)

    for row in iterate.tabular_txt(
            "http://geneontology.org/gene-associations/"
            f"goa_{ORGANISM['files'][organism]}.gaf.gz",
            skiprows=41,
            delimiter="\t",
            usecols=[0, 1, 4, 8, 12]):
        if row[0] == "UniProtKB" and row[3] in namespaces and row[4].split(
                ":")[-1] == str(organism):
            for protein in primary_accession.get(row[1], {row[1]}):
                yield (protein, row[2])

    for row in iterate.tabular_txt(
            "http://geneontology.org/gene-associations/"
            f"goa_{ORGANISM['files'][organism]}_isoform.gaf.gz",
            skiprows=41,
            delimiter="\t",
            usecols=[0, 4, 8, 12, 16]):
        if row[0] == "UniProtKB" and row[2] in namespaces and row[3].split(
                ":")[-1] == str(organism) and row[4].startswith("UniProtKB:"):
            yield (row[4].split(":")[1], row[1])


def convert_namespaces(namespaces: Container[str]) -> tuple[str]:
    """
    Converts Gene Ontology namespace identifiers.

    Args:
        namespaces: Gene Ontology namespaces.

    Returns:
        The corresponding identifiers used in annotation files.
    """
    return tuple({
        "cellular_component": "C",
        "molecular_function": "F",
        "biological_process": "P"
    }[ns] for ns in namespaces)
