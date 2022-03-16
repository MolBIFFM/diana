"""The interface for the Reactome database."""
from typing import Container, Generator, Optional

from access import iterate

from databases import uniprot

ORGANISM = {"data": {9606: "Homo sapiens"}, "file": {9606: "homo_sapiens",}}


def get_protein_protein_interactions(
    interaction_type: Optional[Container[str]] = None,
    interaction_context: Optional[Container[str]] = None,
    taxonomy_identifier: int = 9606,
) -> Generator[tuple[str, str, float], None, None]:
    """
    Yields protein-protein interactions from Reactome.

    Args:
        interaction_type: The accepted interaction type annotation. If none are
            specified, any is accepted.
        interaction_context: The accepted interaction context annotation. If
            none are specified, any is accepted.
        taxonomy_identifier: The taxonomy identifier.

    Yields:
        Pairs of interacting proteins.
    """
    primary_accession = uniprot.get_primary_accession(taxonomy_identifier)

    for row in iterate.tabular_txt(
            "https://reactome.org/download/current/interactors/reactome.{organism}.interactions.tab-delimited.txt"
            .format(organism=ORGANISM["file"].get(taxonomy_identifier,
                                                  "all_species")),
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

                for primary_interactor_a in primary_accession.get(
                        interactor_a, {interactor_a}):
                    for primary_interactor_b in primary_accession.get(
                            interactor_b, {interactor_b}):
                        yield (primary_interactor_a, primary_interactor_b)


def get_pathways(
        taxonomy_identifier: int = 9606
) -> Generator[tuple[str, str], None, None]:
    """
    Yields Reactome pathways.

    Args:
        taxonomy_identifier: The taxonomy identifier.

    Yields:
        Pairs of stable pathway identifier and pathway name.
    """
    for row in iterate.tabular_txt(
            "https://reactome.org/download/current/ReactomePathways.txt",
            delimiter="\t",
            usecols=[0, 1, 2]):
        if not taxonomy_identifier or row[2] == ORGANISM["data"].get(
                taxonomy_identifier):
            yield (row[0], row[1])


def get_pathway_relations() -> Generator[tuple[str, str], None, None]:
    """
    Yields Reactome pathway relations.

    Yields:
        Pairs of parent and child stable pathway identifiers.
    """
    for row in iterate.tabular_txt(
            "https://reactome.org/download/current/ReactomePathwaysRelation.txt",
            delimiter="\t",
            usecols=[0, 1]):
        yield (row[0], row[1])


def get_pathway_map(
        taxonomy_identifier: int = 9606
) -> Generator[tuple[str, str], None, None]:
    """
    Yields Reactome pathway annotations.

    Args:
        taxonomy_identifier: The taxonomy identifier.

    Yields:
        Pairs of protein accession and stable pathway identifier.
    """
    primary_accession = uniprot.get_primary_accession(taxonomy_identifier)

    for row in iterate.tabular_txt(
            "https://reactome.org/download/current/UniProt2Reactome_All_Levels.txt",
            delimiter="\t",
            usecols=[0, 1, 5]):
        if not taxonomy_identifier or row[5] == ORGANISM["data"].get(
                taxonomy_identifier):
            for protein in primary_accession.get(row[0], {row[0]}):
                yield protein, row[1]
