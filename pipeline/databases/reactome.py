"""The interface for the Reactome database."""
from typing import Container, Generator, Optional

from access import iterate

from databases import uniprot

ORGANISM = {"data": {9606: "Homo sapiens"}, "files": {9606: "homo_sapiens"}}


def get_protein_protein_interactions(
    interaction_type: Optional[Container[str]] = None,
    interaction_context: Optional[Container[str]] = None,
    organism: int = 9606,
) -> Generator[tuple[str, str, float], None, None]:
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


def get_pathways(
        organism: Optional[int] = None
) -> Generator[tuple[str, str], None, None]:
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
        if not organism or row[2] == ORGANISM["data"].get(organism):
            yield (row[0], row[1])


def get_pathway_relations() -> Generator[tuple[str, str], None, None]:
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


def get_pathway_annotation(
        organism: Optional[int] = None
) -> Generator[tuple[str, str], None, None]:
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
        if not organism or row[2] == ORGANISM["data"].get(organism):
            for protein in primary_accession.get(row[0], {row[0]}):
                yield protein, row[1]
