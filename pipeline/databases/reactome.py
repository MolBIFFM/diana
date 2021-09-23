from pipeline.utilities import download, uniprot

REACTOME = "https://reactome.org/download/current/interactors/reactome.{organism}.interactions.tab-delimited.txt"

ORGANISM = {
    9606: "homo_sapiens",
}


def add_proteins(
    network,
    interaction_type=[],
    interaction_context=[],
    taxon_identifier=9606,
):
    primary_accession = uniprot.get_primary_accession(taxon_identifier)

    nodes_to_add = set()
    for row in download.tabular_txt(
            REACTOME.format(
                organism=ORGANISM.get(taxon_identifier, "all_species")),
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
                for int_a in primary_accession.get(interactor_a,
                                                   {interactor_a}):
                    for int_b in primary_accession.get(interactor_b,
                                                       {interactor_b}):
                        if (int_a in network and int_b not in network):
                            nodes_to_add.add(int_b)

                        elif (int_a not in network and int_b in network):
                            nodes_to_add.add(int_a)

    network.add_nodes_from(nodes_to_add)


def add_interactions(
    network,
    interaction_type=[],
    interaction_context=[],
    taxon_identifier=9606,
):
    primary_accession = uniprot.get_primary_accession(taxon_identifier,
                                                      network)

    for row in download.tabular_txt(
            REACTOME.format(organism=ORGANISM[taxon_identifier]),
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

                for int_a in primary_accession.get(interactor_a,
                                                   {interactor_a}):
                    for int_b in primary_accession.get(interactor_b,
                                                       {interactor_b}):
                        if int_a in network and int_b in network and int_a != int_b:
                            network.add_edge(int_a, int_b)
                            network.edges[int_a, int_b]["Reactome"] = 0.5
