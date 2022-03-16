"""The interface for the STRING database."""
from typing import Generator
from databases import uniprot
from access import iterate


def get_protein_protein_interactions(
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
        taxonomy_identifier: int = 9606,
        version: float = 11.5) -> Generator[tuple[str, str, float], None, None]:
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
        taxonomy_identifier: The taxonomy identifier.
        version: The version of the STRING database to query.

    Yields:
        Pairs of interacting proteins and the combined STRING score associated
            with the interaction.
    """
    uniprot_id_map = {}
    for row in iterate.tabular_txt(
            "https://stringdb-static.org/download/protein.aliases.v{version}/{taxonomy_identifier}.protein.aliases.v{version}.txt.gz"
            .format(taxonomy_identifier=taxonomy_identifier, version=version),
            delimiter="\t",
            skiprows=1,
            usecols=[0, 1, 2],
    ):
        if "UniProt_AC" in row[2]:
            if row[0] not in uniprot_id_map:
                uniprot_id_map[row[0]] = set()
            uniprot_id_map[row[0]].add(row[1])

    thresholds = {
        column: threshold for column, threshold in {
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

    primary_accession = uniprot.get_primary_accession(taxonomy_identifier)

    for row in iterate.tabular_txt(
            "https://stringdb-static.org/download/protein.physical.links.full.v{version}/{taxonomy_identifier}.protein.links.full.v{version}.txt.gz"
            .format(taxonomy_identifier=taxonomy_identifier, version=version)
            if physical else
            "https://stringdb-static.org/download/protein.links.full.v{version}/{taxonomy_identifier}.protein.links.full.v{version}.txt.gz"
            .format(taxonomy_identifier=taxonomy_identifier, version=version),
            delimiter=" ",
            header=0,
            usecols=["protein1", "protein2"] + list(thresholds.keys()),
    ):
        if all(row[column] / 1000 >= thresholds[column]
               for column in thresholds):
            for interactor_a in uniprot_id_map.get(row["protein1"], set()):
                for interactor_b in uniprot_id_map.get(row["protein2"], set()):
                    for primary_interactor_a in primary_accession.get(
                            interactor_a, {interactor_a}):
                        for primary_interactor_b in primary_accession.get(
                                interactor_b, {interactor_b}):
                            yield (primary_interactor_a, primary_interactor_b,
                                   row["combined_score"] / 1000)
