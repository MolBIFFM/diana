from pipeline.utilities import download, uniprot

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
    primary_accession = uniprot.get_primary_accession()

    nodes_to_add = set()
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
                for a in primary_accession.get(interactor_a, {interactor_a}):
                    for b in primary_accession.get(interactor_b,
                                                   {interactor_b}):
                        if (a in network and b not in network):
                            nodes_to_add.add(b)

                        elif (a not in network and b in network):
                            nodes_to_add.add(a)

    network.add_nodes_from(nodes_to_add)


def add_interactions(
    network,
    interaction_type=[],
    interaction_context=[],
    taxon_identifier=9606,
):
    primary_accession = uniprot.get_primary_accession(network)

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

            if ((not interaction_type
                 or row["Interaction type"] in interaction_type) and
                (not interaction_context
                 or row["Interaction context"] in interaction_context)):

                for a in primary_accession.get(interactor_a, {interactor_a}):
                    for b in primary_accession.get(interactor_b,
                                                   {interactor_b}):
                        if a in network and b in network and a != b:
                            network.add_edge(a, b)
                            network.edges[a, b]["Reactome"] = 0.5
