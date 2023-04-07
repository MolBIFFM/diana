"""The interface for the Reactome database."""
import os
from typing import (Callable, Container, Hashable, Iterable, Iterator, Mapping,
                    Optional, Sequence)

import scipy.stats
from access import iterate
from algorithms import correction
from databases import uniprot

ORGANISM: dict[str, dict[int, str]] = {
    "data": {
        9606: "Homo sapiens"
    },
    "files": {
        9606: "homo_sapiens"
    }
}


def get_protein_interactions(
        interaction_type: Optional[Container[str]] = None,
        interaction_context: Optional[Container[str]] = None,
        organism: int = 9606,
        file: Optional[str] = None,
        file_uniprot: Optional[str] = None) -> Iterator[tuple[str, str]]:
    """
    Yields protein-protein interactions from Reactome.

    Args:
        interaction_type: The accepted interaction type annotation. If none are
            specified, any are accepted.
        interaction_context: The accepted interaction context annotation. If
            none are specified, any are accepted.
        organism: The NCBI taxonomy identifier for the organism of interest.
        file: The optional local file location to parse interactions from.
        file_uniprot: The optional local file location to parse accessions from.

    Yields:
        Pairs of interacting proteins.
    """
    primary_accession = uniprot.get_primary_accession(organism, file_uniprot)

    for row in iterate.tabular_txt(
            "https://reactome.org/download/current/interactors/"
            f"reactome.{ORGANISM['files'][organism]}.interactions."
            "tab-delimited.txt"
            if file is None or not os.path.isfile(file) else file,
            delimiter="\t",
            header=0,
            usecols=[
                "# Interactor 1 uniprot id",
                "Interactor 2 uniprot id",
                "Interaction type",
                "Interaction context",
            ],
    ):
        if (row["# Interactor 1 uniprot id"].split(":")[0] == "uniprotkb" and
                row["Interactor 2 uniprot id"].split(":")[0] == "uniprotkb"):
            interactor_a = row["# Interactor 1 uniprot id"].split(":")[1]
            interactor_b = row["Interactor 2 uniprot id"].split(":")[1]

            if ((not interaction_type or
                 row["Interaction type"] in interaction_type) and
                (not interaction_context or
                 row["Interaction context"] in interaction_context)):

                if "-" in interactor_a and not interactor_a.split(
                        "-")[1].isnumeric():
                    interactor_a = interactor_a.split("-")[0]

                if "-" in interactor_b and not interactor_b.split(
                        "-")[1].isnumeric():
                    interactor_b = interactor_b.split("-")[0]

                for primary_interactor_a in primary_accession.get(
                        interactor_a.split("-")[0], {interactor_a}):
                    for primary_interactor_b in primary_accession.get(
                            interactor_b.split("-")[0], {interactor_b}):
                        yield (primary_interactor_a, primary_interactor_b)


def get_pathways(
        organism: int = 9606,
        file_pathways: Optional[str] = None) -> Iterator[tuple[str, str]]:
    """
    Yields Reactome pathways.

    Args:
        organism: The NCBI taxonomy identifier for the organism of interest.
        file_pathways: The optional local file location to parse pathways
            from.

    Yields:
        Pairs of stable pathway identifier and pathway name.
    """
    for row in iterate.tabular_txt(
            "https://reactome.org/download/current/ReactomePathways.txt"
            if file_pathways is None or not os.path.isfile(file_pathways) else
            file_pathways,
            delimiter="\t",
            usecols=[0, 1, 2]):
        if row[2] == ORGANISM["data"][organism]:
            yield (row[0], row[1])


def get_pathway_relations(
        file_pathways_relation: Optional[str] = None
) -> Iterator[tuple[str, str]]:
    """
    Yields Reactome pathway relations.

    Args:
        file_pathways_relation: The optional local file location to parse
            pathway relations from.

    Yields:
        Pairs of parent and child stable pathway identifiers.
    """
    for row in iterate.tabular_txt(
            "https://reactome.org/download/current/"
            "ReactomePathwaysRelation.txt" if file_pathways_relation is None or
            not os.path.isfile(file_pathways_relation) else
            file_pathways_relation,
            delimiter="\t",
            usecols=[0, 1]):
        yield (row[0], row[1])


def get_pathway_annotation(
        organism: int = 9606,
        file_accession_map: Optional[str] = None,
        file_uniprot: Optional[str] = None) -> Iterator[tuple[str, str]]:
    """
    Yields Reactome pathway annotations.

    Args:
        organism: The NCBI taxonomy identifier for the organism of interest.
        file_accession_map: The optional local file location to parse accession
            associations from.
        file_uniprot: The optional local file location to parse accessions from.

    Yields:
        Pairs of protein accession and stable pathway identifier.
    """
    primary_accession = uniprot.get_primary_accession(organism, file_uniprot)

    for row in iterate.tabular_txt(
            "https://reactome.org/download/current/"
            "UniProt2Reactome_All_Levels.txt" if file_accession_map is None or
            not os.path.isfile(file_accession_map) else file_accession_map,
            delimiter="\t",
            usecols=[0, 1, 5]):
        if row[2] == ORGANISM["data"][organism]:
            for protein in primary_accession.get(row[0], {row[0]}):
                yield protein, row[1]


def get_enrichment(
    proteins: Sequence[frozenset[str]],
    reference: Sequence[Iterable[str]],
    enrichment_test: Callable[[int, int, int, int],
                              float] = lambda k, M, n, N: float(
                                  scipy.stats.hypergeom.sf(k - 1, M, n, N)),
    multiple_testing_correction: Callable[[dict[Hashable, float]], Mapping[
        Hashable, float]] = correction.benjamini_yekutieli,
    organism: int = 9606,
    file_pathways: Optional[str] = None,
    file_accession_map: Optional[str] = None,
    file_uniprot: Optional[str] = None
) -> dict[frozenset[str], dict[tuple[str, str], tuple[float, frozenset[str]]]]:
    """
    Test sets of proteins for enrichment of Reactome pathways.

    Args:
        proteins: The sets of proteins to test.
        reference: Optional reference sets of proteins with respect to which
            enrichment is computed. If not provided, the entire Reactome pathway
            annotation specific to the organism of interest is used.
        enrichment_test: The statistical test used to assess enrichment of
            Reactome pathways.
        multiple_testing_correction: The procedure to correct for multiple
            testing of multiple pathways and sets of proteins.
        organism: The NCBI taxonomy identifier for the organism of interest.
        file_pathways: The optional local file location to parse pathways
            from.
        file_accession_map: The optional local file location to parse accession
            associations from.
        file_uniprot: The optional local file location to parse accessions from.

    Returns:
        Corrected p-value for the enrichment of each Reactome pathway by
        each network and associated proteins.
    """
    name = {}
    for pathway, pathway_name in get_pathways(organism, file_pathways):
        name[pathway] = pathway_name

    annotations: dict[str, set[str]] = {}
    for protein, pathway in get_pathway_annotation(organism, file_accession_map,
                                                   file_uniprot):
        if any(not reference[i if len(reference) == len(proteins) else 0] or
               protein in reference[i if len(reference) == len(proteins) else 0]
               for i in range(len(reference))):
            if pathway not in annotations:
                annotations[pathway] = set()
            annotations[pathway].add(protein)

    annotations = {
        pathway: annotation
        for pathway, annotation in annotations.items()
        if annotations
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
            pathway: annotation.intersection(
                reference[i if len(reference) ==
                          len(proteins) else 0]).intersection(prt)
            if reference[i if len(reference) == len(proteins) else 0] else
            annotation.intersection(prt)
            for pathway, annotation in annotations.items()
        } for i, prt in enumerate(proteins)
    }

    p_value = multiple_testing_correction({(prt, pathway): enrichment_test(
        len(prt_intersection[prt][pathway]),
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
        len(annotated_prt[prt])) for pathway, annotation in annotations.items()
                                           for i, prt in enumerate(proteins)})

    return {
        prt: {(pathway, name[pathway]):
              (p_value[(prt, pathway)],
               frozenset(prt_intersection[prt][pathway]))
              for pathway in annotations} for prt in proteins
    }
