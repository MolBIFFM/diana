"""The interface for the STRING database."""
from typing import Iterator

from access import iterate

from databases import uniprot


def get_protein_interactions(
        neighborhood: float = 0.0,
        neighborhood_transferred: float = 0.0,
        fusion: float = 0.0,
        cooccurence: float = 0.0,
        homology: float = 0.0,
        coexpression: float = 0.0,
        coexpression_transferred: float = 0.0,
        experiments: float = 0.0,
        experiments_transferred: float = 0.0,
        database: float = 0.0,
        database_transferred: float = 0.0,
        textmining: float = 0.0,
        textmining_transferred: float = 0.0,
        combined_score: float = 0.0,
        physical: bool = False,
        organism: int = 9606,
        version: float = 11.5) -> Iterator[tuple[str, str, float]]:
    """
    Yields protein-protein interactions from STRING.

    Args:
        neighborhood: The normal gene neighborhood score threshold.
        neighborhood_transferred: The transferred gene neighborhood score
            threshold.
        fusion: The gene fusion score threshold.
        cooccurrence: The gene cooccurrence score threshold.
        homology: The homology score threshold.
        coexpression: The normal coexpression score threshold.
        coexpression_transferred: The transferred coexpression score threshold.
        experiments: The normal experiments score threshold.
        experiments_transferred: The transferred experiments score threshold.
        database: The normal database score threshold.
        database_transferred: The transferred database score threshold.
        textmining: The normal textmining score threshold.
        textmining_transferred: The transferred textmining score threshold.
        combined_score: The combined score threshold.
        physical: If True, yield only physical interactions.
        organism: The NCBI taxonomy identifier for the organism of interest.
        version: The version of the STRING database.

    Yields:
        Pairs of interacting proteins and the combined STRING score associated
            with the interaction.
    """
    uniprot_id: dict[str, set[str]] = {}
    for row in iterate.tabular_txt(
            f"https://stringdb-static.org/download/protein.aliases.v{version}/"
            f"{organism}.protein.aliases.v{version}.txt.gz",
            delimiter="\t",
            skiprows=1,
            usecols=[0, 1, 2],
    ):
        if "UniProt_AC" in row[2]:
            if row[0] not in uniprot_id:
                uniprot_id[row[0]] = set()
            uniprot_id[row[0]].add(row[1])

    thresholds = {
        column: threshold
        for column, threshold in {
            "neighborhood": neighborhood,
            "neighborhood_transferred": neighborhood_transferred,
            "fusion": fusion,
            "cooccurence": cooccurence,
            "homology": homology,
            "coexpression": coexpression,
            "coexpression_transferred": coexpression_transferred,
            "experiments": experiments,
            "experiments_transferred": experiments_transferred,
            "database": database,
            "database_transferred": database_transferred,
            "textmining": textmining,
            "textmining_transferred": textmining_transferred,
        }.items() if threshold
    }
    thresholds["combined_score"] = combined_score

    primary_accession = uniprot.get_primary_accession(organism)

    for row in iterate.tabular_txt(
            "https://stringdb-static.org/download/"
            f"protein.physical.links.full.v{version}/"
            f"{organism}.protein.links.full.v{version}.txt.gz"
            if physical else "https://stringdb-static.org/download/"
            f"protein.links.full.v{version}/"
            f"{organism}.protein.links.full.v{version}.txt.gz",
            delimiter=" ",
            header=0,
            usecols=["protein1", "protein2"] + list(thresholds),
    ):
        if all(row[column] / 1000 >= thresholds[column]
               for column in thresholds):
            for interactor_a in uniprot_id.get(row["protein1"], set()):
                for interactor_b in uniprot_id.get(row["protein2"], set()):
                    for primary_interactor_a in primary_accession.get(
                            interactor_a, {interactor_a}):
                        for primary_interactor_b in primary_accession.get(
                                interactor_b, {interactor_b}):
                            yield (primary_interactor_a, primary_interactor_b,
                                   row["combined_score"] / 1000)
