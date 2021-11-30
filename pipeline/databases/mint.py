from formats import mitab
from uniprot import uniprot
from download import download

MINT = "http://www.ebi.ac.uk/Tools/webservices/psicquic/mint/webservices/current/search/query/{organism}"

ORGANISM = {9606: "species:human"}


def add_proteins(network,
                 interaction_detection_methods=[],
                 interaction_types=[],
                 mi_score=0.0,
                 taxon_identifier=9606):
    primary_accession = uniprot.get_primary_accession(taxon_identifier)

    nodes_to_add = set()
    for row in download.tabular_txt(
            MINT.format(organism=ORGANISM.get(taxon_identifier, "*")),
            delimiter="\t",
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
        if (int(mitab.get_identifiers_from_namespace(
                row[5], "taxid")) == taxon_identifier
                and int(mitab.get_identifiers_from_namespace(row[6], "taxid"))
                == taxon_identifier and (not interaction_detection_methods
                                         or mitab.namespace_has_any_term_from(
                                             row[6],
                                             "psi-mi",
                                             interaction_detection_methods,
                                         )) and
            (not interaction_types or mitab.namespace_has_any_term_from(
                row[7], "psi-mi", interaction_types))
                and (score := mitab.get_identifiers_from_namespace(
                    row[8], "intact-miscore"))
                and float(score[0]) >= mi_score):
            for interactor_a in mitab.get_identifiers_from_namespace(
                    row[0],
                    "uniprotkb") + mitab.get_identifiers_from_namespace(
                        row[2], "uniprotkb"):
                for interactor_b in mitab.get_identifiers_from_namespace(
                        row[1],
                        "uniprotkb") + mitab.get_identifiers_from_namespace(
                            row[3], "uniprotkb"):
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
            MINT.format(organism=ORGANISM.get(taxon_identifier, "*")),
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
        if ((not interaction_detection_methods
             or mitab.namespace_has_any_term_from(
                 row[4],
                 "psi-mi",
                 interaction_detection_methods,
             )) and
            (not interaction_types or mitab.namespace_has_any_term_from(
                row[7], "psi-mi", interaction_types))
                and (score := mitab.get_identifiers_from_namespace(
                    row[8], "intact-miscore"))
                and float(score[0]) >= mi_score):
            for interactor_a in mitab.get_identifiers_from_namespace(
                    row[0],
                    "uniprotkb") + mitab.get_identifiers_from_namespace(
                        row[2], "uniprotkb"):
                for interactor_b in mitab.get_identifiers_from_namespace(
                        row[1],
                        "uniprotkb") + mitab.get_identifiers_from_namespace(
                            row[3], "uniprotkb"):
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
                                        primary_interactor_b]["MINT"] = max(
                                            float(score[0]),
                                            network.edges[
                                                primary_interactor_a,
                                                primary_interactor_b].get(
                                                    "MINT", 0.0),
                                        )
                                else:
                                    network.add_edge(primary_interactor_a,
                                                     primary_interactor_b)
                                    network.edges[
                                        primary_interactor_a,
                                        primary_interactor_b]["MINT"] = float(
                                            score[0])
