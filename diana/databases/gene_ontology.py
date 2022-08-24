"""The interface for the Gene Ontology database."""
from typing import (Callable, Collection, Container, Hashable, Iterable,
                    Iterator, Mapping)

import scipy.stats
from access import iterate
from algorithms import correction

from databases import uniprot

ORGANISM: dict[str, dict[int, str]] = {"files": {9606: "human"}}


def get_ontology(namespaces: Container[str] = (
    "cellular_component", "molecular_function",
    "biological_process")) -> Iterator[dict[str, str | tuple[str, ...]]]:
    """
    Yields Gene Ontology terms from the given namespaces.

    Args:
        namespaces: The Gene Ontology namespaces to consider terms from.

    Yields:
        Mappings containing a Gene Ontology terms' GO ID, name, namespace,
            related terms and alternative GO IDs.
    """
    term: dict[str, str] = {
        "id": "",
        "name": "",
        "namespace": "",
    }

    is_a: list[str] = []
    alt_id: list[str] = []

    for line in iterate.txt("http://purl.obolibrary.org/obo/go.obo"):
        if line in ("[Term]", "[Typedef]"):
            if term.get("id") and term.get("namespace") in namespaces:
                yield {
                    **term,
                    **{
                        "is_a": tuple(is_a)
                    },
                    **{
                        "alt_id": tuple(alt_id)
                    }
                }

            for attribute in ("id", "name", "namespace"):
                term[attribute] = ""

            is_a.clear()
            alt_id.clear()

        elif any(
                line.startswith(f"{tag}:")
                for tag in ("id", "name", "namespace")):
            term[line.split(":")[0]] = line.split(":", maxsplit=1)[1].strip()

        elif line.startswith("is_a:"):
            is_a.append(line.split(":", maxsplit=1)[1].split("!")[0].strip())

        elif line.startswith("alt_id:"):
            alt_id.append(line.split(":", maxsplit=1)[1].strip())


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


def convert_namespaces(namespaces: Iterable[str]) -> tuple[str, ...]:
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


def get_enrichment(
    proteins: Iterable[frozenset[str]],
    enrichment_test: Callable[
        [int, int, int, int],
        float] = lambda k, M, n, N: scipy.stats.hypergeom.sf(k - 1, M, n, N),
    multiple_testing_correction: Callable[[dict[Hashable, float]], Mapping[
        Hashable, float]] = correction.benjamini_yekutieli,
    organism: int = 9606,
    namespaces: Collection[str] = ("cellular_component", "molecular_function",
                                   "biological_process"),
    annotation_as_reference: bool = False
) -> dict[frozenset, dict[tuple[str, str], float]]:
    """
    Test sets of proteins for enrichment of Gene Ontology terms.

    Args:
        proteins: The sets of proteins test.
        enrichment_test: The statistical test used to assess enrichment of Gene
            Ontology terms.
        multiple_testing_correction: The procedure to correct for multiple
            testing of multiple terms and sets of proteins.
        organism: The NCBI taxonomy identifier for the organism of interest.
        namespaces: The Gene Ontology namespaces.
        annotation_as_reference: If True, compute enrichment with respect to the
            species-specific Gene Ontology annotation in namespaces, otherwise
            with respect to the union of the sets of proteins.

    Returns:
        Corrected p-value for the enrichment of each Gene Ontology term by each
        network.
    """
    name = {}
    go_id: dict[str, set[str]] = {}
    for term in get_ontology(namespaces):
        if isinstance(term["id"], str) and isinstance(term["name"], str):
            name[term["id"]] = term["name"]
            for alt_id in term["alt_id"]:
                if alt_id not in go_id:
                    go_id[alt_id] = set()
                go_id[alt_id].add(term["id"])

    annotation: dict[str, set[str]] = {}
    for protein, annotated_term in get_annotation(
            organism, convert_namespaces(namespaces)):
        if annotation_as_reference or any(protein in prt for prt in proteins):
            for primary_term in go_id.get(annotated_term, {annotated_term}):
                if primary_term not in annotation:
                    annotation[primary_term] = set()
                annotation[primary_term].add(protein)

    annotation = {
        term: proteins for term, proteins in annotation.items() if proteins
    }

    annotated_proteins = set.union(*annotation.values())

    annotated_network_proteins = {
        prt: len(annotated_proteins.intersection(prt)) for prt in proteins
    }

    network_intersection = {
        prt:
        {term: len(annotation[term].intersection(prt)) for term in annotation}
        for prt in proteins
    }

    p_value = multiple_testing_correction({
        (prt, term): enrichment_test(network_intersection[prt][term],
                                     len(annotated_proteins),
                                     len(annotation[term]),
                                     annotated_network_proteins[prt])
        for term in annotation for prt in proteins
    })

    return {
        prt: {(term, name[term]): p_value[(prt, term)] for term in annotation
             } for prt in proteins
    }
