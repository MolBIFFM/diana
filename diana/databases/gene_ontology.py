"""The interface for the Gene Ontology database."""
from typing import (Callable, Collection, Container, Hashable, Iterable,
                    Iterator, Literal, Mapping, Optional, Sequence)

import scipy.stats
from access import iterate
from algorithms import correction
from databases import uniprot

ORGANISM: dict[str, dict[int, str]] = {"files": {9606: "human"}}


def get_ontology(
    namespaces: Container[Literal["cellular_component", "molecular_function",
                                  "biological_process"]] = (),
    file_ontology: Optional[str] = None
) -> Iterator[dict[str, str | tuple[str, ...]]]:
    """
    Yields Gene Ontology terms from the given namespaces.

    Args:
        namespaces: The Gene Ontology namespaces to consider terms from. If
            empty, any namespace is considered.
        file_ontology: The optional local file location to parse terms from.

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

    for line in iterate.txt("http://purl.obolibrary.org/obo/go.obo"
                            if file_ontology is None else file_ontology):
        if line in ("[Term]", "[Typedef]"):
            if (not namespaces or
                    term.get("id") and term.get("namespace") in namespaces):
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
        namespaces: Container[Literal["C", "F", "P"]] = (),
        file_annotation: Optional[str] = None,
        file_annotation_isoform: Optional[str] = None,
        file_uniprot: Optional[str] = None) -> Iterator[tuple[str, str]]:
    """
    Yields Gene Ontology annotations within specified namespaces.

    Args:
        organism: The NCBI taxonomy identifier for the organism of interest.
        namespace: The Gene Ontology namespace identifiers. If empty, any
            namespace is considered.
        file_annotation: The optional local file location to parse annotations
            from.
        file_annotation: The optional local file location to parse isoform
            annotations from.
        file_uniprot: The optional local file location to parse accessions from.

    Yields:
        Pairs of protein accessions and Gene Ontology term identifiers.
    """
    primary_accession = uniprot.get_primary_accession(organism, file_uniprot)

    for row in iterate.tabular_txt(
            "http://geneontology.org/gene-associations/"
            f"goa_{ORGANISM['files'][organism]}.gaf.gz"
            if file_annotation is None else file_annotation,
            skiprows=41,
            delimiter="\t",
            usecols=[0, 1, 4, 8, 12]):
        if (row[0] == "UniProtKB" and
            (not namespaces or row[3] in namespaces) and
                row[4].split(":")[-1] == str(organism)):
            for protein in primary_accession.get(row[1], {row[1]}):
                yield (protein, row[2])

    for row in iterate.tabular_txt(
            "http://geneontology.org/gene-associations/"
            f"goa_{ORGANISM['files'][organism]}_isoform.gaf.gz"
            if file_annotation_isoform is None else file_annotation_isoform,
            skiprows=41,
            delimiter="\t",
            usecols=[0, 4, 8, 12, 16]):
        if (row[0] == "UniProtKB" and
            (not namespaces or row[2] in namespaces) and
                row[3].split(":")[-1] == str(organism) and
                row[4].startswith("UniProtKB:")):
            yield (row[4].split(":")[1], row[1])


def convert_namespaces(
    namespaces: Iterable[Literal["cellular_component", "molecular_function",
                                 "biological_process"]]
) -> tuple[Literal["C", "F", "P"], ...]:
    """
    Converts Gene Ontology namespace identifiers.

    Args:
        namespaces: Gene Ontology namespaces.

    Returns:
        The corresponding identifiers used in annotation files.
    """

    conversion: dict[Literal["cellular_component", "molecular_function",
                             "biological_process"], Literal["C", "F", "P"]] = {
                                 "cellular_component": "C",
                                 "molecular_function": "F",
                                 "biological_process": "P"
                             }

    return tuple(conversion[ns] for ns in namespaces)


def get_enrichment(
    proteins: Sequence[frozenset[str]],
    reference: Sequence[Iterable[str]],
    enrichment_test: Callable[[int, int, int, int],
                              float] = lambda k, M, n, N: float(
                                  scipy.stats.hypergeom.sf(k - 1, M, n, N)),
    multiple_testing_correction: Callable[[dict[Hashable, float]], Mapping[
        Hashable, float]] = correction.benjamini_yekutieli,
    organism: int = 9606,
    namespaces: Collection[Literal["cellular_component", "molecular_function",
                                   "biological_process"]] = frozenset(),
    file_ontology: Optional[str] = None,
    file_annotation: Optional[str] = None,
    file_annotation_isoform: Optional[str] = None,
    file_uniprot: Optional[str] = None
) -> dict[frozenset[str], dict[tuple[str, str], tuple[float, frozenset[str]]]]:
    """
    Test sets of proteins for enrichment of Gene Ontology terms.

    Args:
        proteins: The sets of proteins test.
        reference: Optional reference sets of proteins with respect to which
            enrichment is computed. If not provided, the entire Gene Ontology
            annotation specific to the organism of interest is used.
        enrichment_test: The statistical test used to assess enrichment of Gene
            Ontology terms.
        multiple_testing_correction: The procedure to correct for multiple
            testing of multiple terms and sets of proteins.
        organism: The NCBI taxonomy identifier for the organism of interest.
        namespaces: The Gene Ontology namespaces to consider. If empty, any
            namespace is considered.
        file_ontology: The optional local file location to parse terms from.
        file_annotation: The optional local file location to parse annotations
            from.
        file_annotation: The optional local file location to parse isoform
            annotations from.
        file_uniprot: The optional local file location to parse accessions from.

    Returns:
        Corrected p-value for the enrichment of each Gene Ontology term by each
        network and associated proteins.
    """
    name = {}
    go_id: dict[str, set[str]] = {}
    for term in get_ontology(namespaces, file_ontology):
        if isinstance(term["id"], str) and isinstance(term["name"], str):
            name[term["id"]] = term["name"]
            for alt_id in term["alt_id"]:
                if alt_id not in go_id:
                    go_id[alt_id] = set()
                go_id[alt_id].add(term["id"])

    annotations: dict[str, set[str]] = {}
    for protein, annotated_term in get_annotation(
            organism, convert_namespaces(namespaces), file_annotation,
            file_annotation_isoform, file_uniprot):
        if any(not reference[i if len(reference) == len(proteins) else 0] or
               protein in reference[i if len(reference) == len(proteins) else 0]
               for i in range(len(reference))):
            for primary_term in go_id.get(annotated_term, {annotated_term}):
                if primary_term not in annotations:
                    annotations[primary_term] = set()
                annotations[primary_term].add(protein)

    annotations = {
        term: annotation
        for term, annotation in annotations.items()
        if annotation
    }

    annotated_proteins = set.union(*annotations.values())

    annotated_prt = {
        prt: annotated_proteins.intersection(
            reference[i if len(reference) ==
                      len(proteins) else 0]).intersection(prt)
        if reference[i if len(reference) == len(proteins) else 0] else
        annotated_proteins.intersection(prt) for i, prt in enumerate(proteins)
    }

    prt_intersection = {
        prt: {
            term: annotation.intersection(
                reference[i if len(reference) ==
                          len(proteins) else 0]).intersection(prt)
            if reference[i if len(reference) == len(proteins) else 0] else
            annotation.intersection(prt)
            for term, annotation in annotations.items()
        } for i, prt in enumerate(proteins)
    }

    p_value = multiple_testing_correction({(prt, term): enrichment_test(
        len(prt_intersection[prt][term]),
        len(
            annotated_proteins.intersection(reference[i if len(reference) ==
                                                      len(proteins) else 0])
            if reference[i if len(reference) ==
                         len(proteins) else 0] else annotated_proteins),
        len(
            annotation.intersection(reference[i if len(reference) ==
                                              len(proteins) else 0])
            if reference[i if len(reference) ==
                         len(proteins) else 0] else annotation),
        len(annotated_prt[prt])) for term, annotation in annotations.items()
                                           for i, prt in enumerate(proteins)})

    return {
        prt: {(term, name[term]):
              (p_value[(prt, term)], frozenset(prt_intersection[prt][term]))
              for term in annotations} for prt in proteins
    }
