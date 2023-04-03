"""The interface for the IntAct database."""
import re
from typing import Container, Iterator, Optional

from access import iterate
from formats import mitab

from databases import uniprot


def get_protein_interactions(
        interaction_detection_methods: Optional[Container[str]] = None,
        interaction_types: Optional[Container[str]] = None,
        psi_mi_score: float = 0.0,
        organism: int = 9606,
        file: Optional[str] = None,
        file_uniprot: Optional[str] = None) -> Iterator[tuple[str, str, float]]:
    """
    Yields protein-protein interactions from IntAct.

    Args:
        interaction_detection_methods: The accepted PSI-MI identifiers or terms
            for interaction detection method. If none are specified, any are
            accepted.
        interaction_types: The accepted PSI-MI identifiers or terms for
            interaction type. If none are specified, any are accepted.
        psi_mi_score: The PSI-MI score threshold.
        organism: The NCBI taxonomy identifier for the organism of interest.
        file: The optional local file location to parse interactions from.
        file_uniprot: The optional local file location to parse accessions from.

    Yields:
        Pairs of interacting proteins and the PSI-MI score associated with the
        interaction.
    """
    primary_accession = uniprot.get_primary_accession(organism, file_uniprot)

    for row in iterate.tabular_txt(
            "ftp://ftp.ebi.ac.uk/pub/databases/intact/current/psimitab/"
            "intact.zip" if file is None else file,
            file_from_zip_archive=re.compile(r"^intact\.txt$"),
            delimiter="\t",
            header=0,
            usecols=[
                "#ID(s) interactor A", "ID(s) interactor B",
                "Alt. ID(s) interactor A", "Alt. ID(s) interactor B",
                "Interaction detection method(s)", "Taxid interactor A",
                "Taxid interactor B", "Interaction type(s)",
                "Confidence value(s)", "Type(s) interactor A",
                "Type(s) interactor B"
            ],
    ):
        if ((not interaction_detection_methods or
             mitab.namespace_has_any_identifier_from(
                 row["Interaction detection method(s)"],
                 "psi-mi",
                 interaction_detection_methods,
             ) or mitab.namespace_has_any_term_from(
                 row["Interaction detection method(s)"],
                 "psi-mi",
                 interaction_detection_methods,
             )) and (mitab.namespace_has_identifier(row["Taxid interactor A"],
                                                    "taxid", organism) and
                     mitab.namespace_has_identifier(row["Taxid interactor B"],
                                                    "taxid", organism)) and
            (not interaction_types or mitab.namespace_has_any_identifier_from(
                row["Interaction type(s)"], "psi-mi", interaction_types) or
             mitab.namespace_has_any_term_from(row["Interaction type(s)"],
                                               "psi-mi", interaction_types)) and
            (score := mitab.get_identifiers_from_namespace(
                row["Confidence value(s)"], "intact-miscore")) and
                float(score[0]) >= psi_mi_score and mitab.namespace_has_term(
                    row["Type(s) interactor A"], "psi-mi", "protein") and
                mitab.namespace_has_term(row["Type(s) interactor B"], "psi-mi",
                                         "protein")):
            for interactor_a in mitab.get_identifiers_from_namespace(
                    row["#ID(s) interactor A"],
                    "uniprotkb") + mitab.get_identifiers_from_namespace(
                        row["Alt. ID(s) interactor A"], "uniprotkb"):

                if "-" in interactor_a and not interactor_a.split(
                        "-")[1].isnumeric():
                    interactor_a = interactor_a.split("-")[0]

                for interactor_b in mitab.get_identifiers_from_namespace(
                        row["ID(s) interactor B"],
                        "uniprotkb") + mitab.get_identifiers_from_namespace(
                            row["Alt. ID(s) interactor B"], "uniprotkb"):

                    if "-" in interactor_b and not interactor_b.split(
                            "-")[1].isnumeric():
                        interactor_b = interactor_b.split("-")[0]

                    for primary_interactor_a in primary_accession.get(
                            interactor_a.split("-")[0], {interactor_a}):
                        for primary_interactor_b in primary_accession.get(
                                interactor_b.split("-")[0], {interactor_b}):
                            yield (primary_interactor_a, primary_interactor_b,
                                   float(score[0]))
