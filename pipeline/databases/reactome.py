import networkx as nx

from pipeline.utilities import download

REACTOME = "https://reactome.org/download/current/interactors/reactome.{organism}.interactions.tab-delimited.txt"

ORGANISM = {
    9606: {
        REACTOME: "homo_sapiens",
    },
}


def add_proteins(
    network,
    interaction_type=[],
    interaction_context=[],
    taxon_identifier=9606,
):
    for row in download.tabular_txt(
            REACTOME.format(organism=ORGANISM[taxon_identifier][REACTOME]),
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

            if (not interaction_type
                    or row["Interaction type"] in interaction_type) and (
                        not interaction_context
                        or row["Interaction context"] in interaction_context):
                if (interactor_a in network and interactor_b not in network
                        and network.nodes[interactor_a].get("protein")):
                    network.add_node(interactor_b)

                elif (interactor_a not in network and interactor_b in network
                      and network.nodes[interactor_b].get("protein")):
                    network.add_node(interactor_a)


def add_interactions(
    network,
    interaction_type=[],
    interaction_context=[],
    taxon_identifier=9606,
):

    for row in download.tabular_txt(
            REACTOME.format(organism=ORGANISM[taxon_identifier][REACTOME]),
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

            if (interactor_a != interactor_b
                    and (not interaction_type
                         or row["Interaction type"] in interaction_type)
                    and (not interaction_context
                         or row["Interaction context"] in interaction_context)
                    and interactor_a in network and interactor_b in network):

                network.add_edge(interactor_a, interactor_b)
                network.edges[interactor_a, interactor_b]["Reactome"] = 0.5
