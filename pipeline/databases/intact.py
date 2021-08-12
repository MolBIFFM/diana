import networkx as nx

from pipeline.utilities import download, mitab

INTACT_ZIP_ARCHIVE = (
    "ftp://ftp.ebi.ac.uk/pub/databases/intact/current/psimitab/intact.zip")
INTACT = "intact.txt"


def add_proteins(
    network,
    interaction_detection_methods=[],
    interaction_types=[],
    mi_score=0.0,
):
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
        nodes_to_add = set()
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
            if (interactor_a in network and interactor_b not in network):
                nodes_to_add.add(interactor_b)

            elif (interactor_a not in network and interactor_b in network):
                nodes_to_add.add(interactor_a)
    network.add_nodes_from(nodes_to_add)


def add_interactions(
    network,
    interaction_detection_methods=[],
    interaction_types=[],
    mi_score=0.0,
):

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
                    row["Confidence value(s)"], "intact-miscore"))
                and interactor_a in network and interactor_b in network):
            score = float(score)

            if score >= mi_score:
                if network.has_edge(interactor_a, interactor_b):
                    network.edges[interactor_a, interactor_b]["IntAct"] = max(
                        score,
                        network.edges[interactor_a,
                                      interactor_b].get("IntAct", 0.0),
                    )
                else:
                    network.add_edge(interactor_a, interactor_b)
                    network.edges[interactor_a, interactor_b]["IntAct"] = score
