from os import PRIO_PGRP
from pipeline.utilities import download, mitab, uniprot

INTACT_ZIP_ARCHIVE = (
    "ftp://ftp.ebi.ac.uk/pub/databases/intact/current/psimitab/intact.zip")
INTACT = "intact.txt"


def add_proteins(
    network,
    interaction_detection_methods=[],
    interaction_types=[],
    mi_score=0.0,
):
    primary_accession = uniprot.get_primary_accession()

    nodes_to_add = set()
    for row in download.tabular_txt(
            INTACT_ZIP_ARCHIVE,
            file_from_zip_archive=INTACT,
            delimiter="\t",
            header=0,
            usecols=[
                "#ID(s) interactor A",
                "ID(s) interactor B",
                "Alt. ID(s) interactor A",
                "Alt. ID(s) interactor B",
                "Taxid interactor A",
                "Taxid interactor B",
                "Interaction detection method(s)",
                "Interaction type(s)",
                "Confidence value(s)",
            ],
    ):
        if ((interactor_a := mitab.get_identifier_from_namespace(
                row["#ID(s) interactor A"], "uniprotkb"))
                and (interactor_b := mitab.get_identifier_from_namespace(
                    row["ID(s) interactor B"], "uniprotkb"))
                and mitab.get_identifier_from_namespace(
                    row["Taxid interactor A"],
                    "taxid") == mitab.get_identifier_from_namespace(
                        row["Taxid interactor B"], "taxid")
                and (not interaction_detection_methods
                     or mitab.namespace_has_any_term_from(
                         row["Interaction detection method(s)"],
                         "psi-mi",
                         interaction_detection_methods,
                     )) and
            (not interaction_types or mitab.namespace_has_any_term_from(
                row["Interaction type(s)"], "psi-mi", interaction_types))
                and (score := mitab.get_identifier_from_namespace(
                    row["Confidence value(s)"], "intact-miscore"))
                and float(score) >= mi_score):
            for a in primary_accession.get(interactor_a, {interactor_a}):
                for b in primary_accession.get(interactor_b, {interactor_b}):
                    if (a in network and b not in network):
                        nodes_to_add.add(b)

                    elif (a not in network and b in network):
                        nodes_to_add.add(a)

    network.add_nodes_from(nodes_to_add)


def add_interactions(
    network,
    interaction_detection_methods=[],
    interaction_types=[],
    mi_score=0.0,
):
    primary_accession = uniprot.get_primary_accession(network)

    for row in download.tabular_txt(
            INTACT_ZIP_ARCHIVE,
            file_from_zip_archive=INTACT,
            delimiter="\t",
            header=0,
            usecols=[
                "#ID(s) interactor A",
                "ID(s) interactor B",
                "Alt. ID(s) interactor A",
                "Alt. ID(s) interactor B",
                "Taxid interactor A",
                "Taxid interactor B",
                "Interaction detection method(s)",
                "Interaction type(s)",
                "Confidence value(s)",
            ],
    ):

        if ((interactor_a := mitab.get_identifier_from_namespace(
                row["#ID(s) interactor A"], "uniprotkb"))
                and (interactor_b := mitab.get_identifier_from_namespace(
                    row["ID(s) interactor B"], "uniprotkb"))
                and interactor_a != interactor_b
                and (not interaction_detection_methods
                     or mitab.namespace_has_any_term_from(
                         row["Interaction detection method(s)"],
                         "psi-mi",
                         interaction_detection_methods,
                     )) and
            (not interaction_types or mitab.namespace_has_any_term_from(
                row["Interaction type(s)"], "psi-mi", interaction_types))
                and (score := mitab.get_identifier_from_namespace(
                    row["Confidence value(s)"], "intact-miscore"))):
            score = float(score)
            if score >= mi_score:
                for a in primary_accession.get(interactor_a, {interactor_a}):
                    for b in primary_accession.get(interactor_b,
                                                   {interactor_b}):
                        if a in network and b in network and a != b:
                            if network.has_edge(a, b):
                                network.edges[a, b]["IntAct"] = max(
                                    score,
                                    network.edges[a, b].get("IntAct", 0.0),
                                )
                            else:
                                network.add_edge(a, b)
                                network.edges[a, b]["IntAct"] = score
