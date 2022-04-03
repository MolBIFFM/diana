"""The interface for the MINT database."""
from typing import Container, Generator, Optional

from access import iterate
from formats import mitab

from databases import uniprot

ORGANISM = {"file": {9606: "species:human"}}


def get_protein_protein_interactions(
    interaction_detection_methods: Optional[Container[str]] = None,
    interaction_types: Optional[Container[str]] = None,
    psi_mi_score: float = 0.0,
    taxonomy_identifier: int = 9606
) -> Generator[tuple[str, str, float], None, None]:
    """
    Yields protein-protein interactions from MINT.

    Args:
        interaction_detection_methods: The accepted PSI-MI terms for interaction
            detection method. If none are specified, any is accepted.
        interaction_types: The accepted PSI-MI terms for interaction type.
            If none are specified, any is accepted.
        psi_mi_score: The PSI-MI score threshold.
        taxonomy_identifier: The taxonomy identifier.

    Yields:
        Pairs of interacting proteins and the PSI-MI score associated with the
        interaction.
    """
    primary_accession = uniprot.get_primary_accession(taxonomy_identifier)

    for row in iterate.tabular_txt(
            f"https://www.ebi.ac.uk/Tools/webservices/psicquic/mint/webservices/current/search/query/{ORGANISM['file'][taxonomy_identifier]}",
            delimiter="\t",
            header=0,
            usecols=[
                0,
                1,
                2,
                3,
                6,
                9,
                10,
                11,
                14,
            ],
    ):
        if ((not interaction_detection_methods or
             mitab.namespace_has_any_term_from(
                 row[4],
                 "psi-mi",
                 interaction_detection_methods,
             )) and (mitab.namespace_has_identifier(row[5], "taxid",
                                                    taxonomy_identifier) and
                     mitab.namespace_has_identifier(row[6], "taxid",
                                                    taxonomy_identifier)) and
            (not interaction_types or mitab.namespace_has_any_term_from(
                row[7], "psi-mi", interaction_types)) and
            (score := mitab.get_identifiers_from_namespace(
                row[8], "intact-miscore")) and float(score[0]) >= psi_mi_score):
            for interactor_a in mitab.get_identifiers_from_namespace(
                    row[0], "uniprotkb") + mitab.get_identifiers_from_namespace(
                        row[2], "uniprotkb"):

                if "-" in interactor_a and not interactor_a.split(
                        "-")[1].isnumeric():
                    interactor_a = interactor_a.split("-")[0]

                for interactor_b in mitab.get_identifiers_from_namespace(
                        row[1],
                        "uniprotkb") + mitab.get_identifiers_from_namespace(
                            row[3], "uniprotkb"):

                    if "-" in interactor_b and not interactor_b.split(
                            "-")[1].isnumeric():
                        interactor_b = interactor_b.split("-")[0]

                    for primary_interactor_a in primary_accession.get(
                            interactor_a, {interactor_a}):
                        for primary_interactor_b in primary_accession.get(
                                interactor_b, {interactor_b}):
                            yield (primary_interactor_a, primary_interactor_b,
                                   float(score[0]))
