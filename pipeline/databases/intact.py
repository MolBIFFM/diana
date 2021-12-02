from formats import mitab
from uniprot import uniprot
from download import download

INTACT_ZIP_ARCHIVE = "ftp://ftp.ebi.ac.uk/pub/databases/intact/current/psimitab/intact.zip"
INTACT = "intact.txt"


def add_proteins(network,
                 interaction_detection_methods=[],
                 interaction_types=[],
                 mi_score=0.0,
                 taxon_identifier=9606):
    primary_accession = uniprot.get_primary_accession(taxon_identifier)

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
                "Interaction type(s)",
                "Interaction detection method(s)",
                "Confidence value(s)",
            ],
    ):
        if (int(
                mitab.get_identifiers_from_namespace(
                    row["Taxid interactor A"], "taxid")) == taxon_identifier
                and int(
                    mitab.get_identifiers_from_namespace(
                        row["Taxid interactor B"],
                        "taxid")) == taxon_identifier
                and (not interaction_detection_methods
                     or mitab.namespace_has_any_term_from(
                         row["Interaction detection method(s)"],
                         "psi-mi",
                         interaction_detection_methods,
                     )) and
            (not interaction_types or mitab.namespace_has_any_term_from(
                row["Interaction type(s)"], "psi-mi", interaction_types))
                and (score := mitab.get_identifiers_from_namespace(
                    row["Confidence value(s)"], "intact-miscore"))
                and float(score[0]) >= mi_score):
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
                            interactor_a, {interactor_a}):
                        for primary_interactor_b in primary_accession.get(
                                interactor_b, {interactor_b}):
                            if (primary_interactor_a in network
                                    and primary_interactor_b not in network):
                                nodes_to_add.add(primary_interactor_b)

                            elif (primary_interactor_a not in network
                                  and primary_interactor_b in network):
                                nodes_to_add.add(primary_interactor_a)

    network.add_nodes_from(nodes_to_add)


def add_interactions(network,
                     interaction_detection_methods=[],
                     interaction_types=[],
                     mi_score=0.0,
                     taxon_identifier=9606):
    primary_accession = uniprot.get_primary_accession(taxon_identifier,
                                                      network)

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
                "Interaction type(s)",
                "Interaction detection method(s)",
                "Confidence value(s)",
            ],
    ):
        if ((not interaction_detection_methods
             or mitab.namespace_has_any_term_from(
                 row["Interaction detection method(s)"],
                 "psi-mi",
                 interaction_detection_methods,
             )) and
            (not interaction_types or mitab.namespace_has_any_term_from(
                row["Interaction type(s)"], "psi-mi", interaction_types))
                and (score := mitab.get_identifiers_from_namespace(
                    row["Confidence value(s)"], "intact-miscore"))
                and float(score[0]) >= mi_score):
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
                            interactor_a, {interactor_a}):
                        for primary_interactor_b in primary_accession.get(
                                interactor_b, {interactor_b}):
                            if (primary_interactor_a in network
                                    and primary_interactor_b in network
                                    and primary_interactor_a !=
                                    primary_interactor_b):
                                if network.has_edge(primary_interactor_a,
                                                    primary_interactor_b):
                                    network.edges[
                                        primary_interactor_a,
                                        primary_interactor_b]["IntAct"] = max(
                                            float(score[0]),
                                            network.edges[
                                                primary_interactor_a,
                                                primary_interactor_b].get(
                                                    "IntAct", 0.0),
                                        )
                                else:
                                    network.add_edge(primary_interactor_a,
                                                     primary_interactor_b)
                                    network.edges[primary_interactor_a,
                                                  primary_interactor_b][
                                                      "IntAct"] = float(
                                                          score[0])
