"""The interface for the Reactome database."""
from typing import Callable, Container, Iterable, Iterator, Mapping, Optional

import scipy.stats
from access import iterate
from analysis import correction

from databases import uniprot

ORGANISM = {"data": {9606: "Homo sapiens"}, "files": {9606: "homo_sapiens"}}


def get_protein_protein_interactions(
    interaction_type: Optional[Container[str]] = None,
    interaction_context: Optional[Container[str]] = None,
    organism: int = 9606,
) -> Iterator[tuple[str, str, float]]:
    """
    Yields protein-protein interactions from Reactome.

    Args:
        interaction_type: The accepted interaction type annotation. If none are
            specified, any are accepted.
        interaction_context: The accepted interaction context annotation. If
            none are specified, any are accepted.
        organism: The NCBI taxonomy identifier for the organism of interest.

    Yields:
        Pairs of interacting proteins.
    """
    primary_accession = uniprot.get_primary_accession(organism)

    for row in iterate.tabular_txt(
            "https://reactome.org/download/current/interactors/"
            f"reactome.{ORGANISM['files'][organism]}.interactions."
            "tab-delimited.txt",
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


def get_pathways(organism: int = 9606) -> Iterator[tuple[str, str]]:
    """
    Yields Reactome pathways.

    Args:
        organism: The NCBI taxonomy identifier for the organism of interest.

    Yields:
        Pairs of stable pathway identifier and pathway name.
    """
    for row in iterate.tabular_txt(
            "https://reactome.org/download/current/ReactomePathways.txt",
            delimiter="\t",
            usecols=[0, 1, 2]):
        if row[2] == ORGANISM["data"][organism]:
            yield (row[0], row[1])


def get_pathway_relations() -> Iterator[tuple[str, str]]:
    """
    Yields Reactome pathway relations.

    Yields:
        Pairs of parent and child stable pathway identifiers.
    """
    for row in iterate.tabular_txt(
            "https://reactome.org/download/current/"
            "ReactomePathwaysRelation.txt",
            delimiter="\t",
            usecols=[0, 1]):
        yield (row[0], row[1])


def get_pathway_annotation(organism: int = 9606) -> Iterator[tuple[str, str]]:
    """
    Yields Reactome pathway annotations.

    Args:
        organism: The NCBI taxonomy identifier for the organism of interest.

    Yields:
        Pairs of protein accession and stable pathway identifier.
    """
    primary_accession = uniprot.get_primary_accession(organism)

    for row in iterate.tabular_txt(
            "https://reactome.org/download/current/"
            "UniProt2Reactome_All_Levels.txt",
            delimiter="\t",
            usecols=[0, 1, 5]):
        if row[2] == ORGANISM["data"][organism]:
            for protein in primary_accession.get(row[0], {row[0]}):
                yield protein, row[1]


def get_enrichment(
    proteins: Iterable[frozenset[str]],
    enrichment_test: Callable[
        [int, int, int, int],
        float] = lambda k, M, n, N: scipy.stats.hypergeom.sf(k - 1, M, n, N),
    multiple_testing_correction: Callable[
        [Mapping[tuple[frozenset[str], str],
                 float]], Mapping[tuple[frozenset[str], str],
                                  float]] = correction.benjamini_hochberg,
    organism: int = 9606,
    annotation_as_reference: bool = True
) -> dict[frozenset[str], dict[tuple[str, str], float]]:
    """
    Test sets of proteins for enrichment of Reactome pathways.

    Args:
        networks: The sets of proteins to test.
        enrichment_test: The statistical test used to assess enrichment of
            Reactome pathways.
        multiple_testing_correction: The procedure to correct for multiple
            testing of multiple pathways and sets of proteins.
        organism: The NCBI taxonomy identifier for the organism of interest.
        annotation_as_reference: If True, compute enrichment with respect to the
            species-specific Reactome pathway annotation, else with respect to
            the union of the sets of proteins.

    Returns:
        Corrected p-value for the enrichment of each Reactome pathway by
        each network.
    """
    name = {}
    for pathway, pathway_name in get_pathways(organism):
        name[pathway] = pathway_name

    annotation = {}
    for protein, pathway in get_pathway_annotation(organism):
        if annotation_as_reference or any(protein in prt for prt in proteins):
            if pathway not in annotation:
                annotation[pathway] = set()
            annotation[pathway].add(protein)

    annotation = {
        pathway: proteins
        for pathway, proteins in annotation.items()
        if proteins
    }

    annotated_proteins = set.union(*annotation.values())

    annotated_network_proteins = {
        prt: len(annotated_proteins.intersection(prt)) for prt in proteins
    }

    network_intersection = {
        prt: {
            pathway: len(annotation[pathway].intersection(prt))
            for pathway in annotation
        } for prt in proteins
    }

    p_value = multiple_testing_correction({
        (prt, pathway): enrichment_test(network_intersection[prt][pathway],
                                        len(annotated_proteins),
                                        len(annotation[pathway]),
                                        annotated_network_proteins[prt])
        for pathway in annotation for prt in proteins
    })

    return {
        prt: {(pathway, name[pathway]): p_value[(prt, pathway)]
              for pathway in annotation} for prt in proteins
    }
