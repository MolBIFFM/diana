from formats import mitab
from uniprot import uniprot
from download import download

DIP = "https://dip.doe-mbi.ucla.edu/dip/File.cgi?FN=2017/tab25/{organism}20170205.txt.gz"

ORGANISM = {9606: "Hsapi"}


def add_proteins(network,
                 interaction_detection_methods=[],
                 interaction_types=[],
                 taxon_identifier=9606):
    primary_accession = uniprot.get_primary_accession(taxon_identifier)

    nodes_to_add = set()
    for row in download.tabular_txt(
            DIP.format(organism=ORGANISM[taxon_identifier]),
            delimiter="\t",
            header=0,
            usecols=[
                "ID interactor A",
                "ID interactor B",
                "Alt. ID interactor A",
                "Alt. ID interactor B",
                "Taxid interactor A",
                "Taxid interactor B",
                "Interaction type(s)",
                "Interaction detection method(s)",
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
                row["Interaction type(s)"], "psi-mi", interaction_types))):
            for interactor_a in mitab.get_identifiers_from_namespace(
                    row["ID interactor A"],
                    "uniprotkb") + mitab.get_identifiers_from_namespace(
                        row["Alt. ID interactor A"], "uniprotkb"):

                if "-" in interactor_a and not interactor_a.split(
                        "-")[1].isnumeric():
                    interactor_a = interactor_a.split("-")[0]

                for interactor_b in mitab.get_identifiers_from_namespace(
                        row["ID interactor B"],
                        "uniprotkb") + mitab.get_identifiers_from_namespace(
                            row["Alt. ID interactor B"], "uniprotkb"):

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
                     taxon_identifier=9606):
    primary_accession = uniprot.get_primary_accession(taxon_identifier,
                                                      network)

    for row in download.tabular_txt(
            DIP.format(organism=ORGANISM[taxon_identifier]),
            delimiter="\t",
            header=0,
            usecols=[
                "ID interactor A",
                "ID interactor B",
                "Alt. ID interactor A",
                "Alt. ID interactor B",
                "Taxid interactor A",
                "Taxid interactor B",
                "Interaction type(s)",
                "Interaction detection method(s)",
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
                row["Interaction type(s)"], "psi-mi", interaction_types))):
            for interactor_a in mitab.get_identifiers_from_namespace(
                    row["ID interactor A"],
                    "uniprotkb") + mitab.get_identifiers_from_namespace(
                        row["Alt. ID interactor A"], "uniprotkb"):

                if "-" in interactor_a and not interactor_a.split(
                        "-")[1].isnumeric():
                    interactor_a = interactor_a.split("-")[0]

                for interactor_b in mitab.get_identifiers_from_namespace(
                        row["ID interactor B"],
                        "uniprotkb") + mitab.get_identifiers_from_namespace(
                            row["Alt. ID interactor B"], "uniprotkb"):

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
                                network.add_edge(primary_interactor_a,
                                                 primary_interactor_b)
                                network.edges[
                                    primary_interactor_a,
                                    primary_interactor_b]["DIP"] = 1.0
