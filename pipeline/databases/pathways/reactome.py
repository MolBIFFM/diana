"""The interface for the Reactome database."""
from typing import Generator

from access import iterate

from databases import uniprot

ORGANISM = {"data": {9606: "Homo sapiens"}}


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
        if not taxonomy_identifier or row[2] == ORGANISM["data"].get(
                taxonomy_identifier):
            for protein in primary_accession.get(row[0], {row[0]}):
                yield protein, row[1]
