from databases import uniprot
from download import download

REACTOME_INTERACTIONS = "https://reactome.org/download/current/interactors/reactome.{organism}.interactions.tab-delimited.txt"
REACTOME_PATHWAYS = "https://reactome.org/download/current/ReactomePathways.txt"
REACTOME_PATHWAY_RELATIONS = "https://reactome.org/download/current/ReactomePathwaysRelation.txt"

ORGANISM = {
    "files": {
        9606: "homo_sapiens",
    },
    "data": {
        9606: "Homo sapiens"
    }
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
            REACTOME_INTERACTIONS.format(organism=ORGANISM["file"].get(
                taxon_identifier, "all_species")),
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


def add_protein_protein_interactions(
    network,
    interaction_type=[],
    interaction_context=[],
    taxon_identifier=9606,
):
    primary_accession = uniprot.get_primary_accession(taxon_identifier,
                                                      network)

    for row in download.tabular_txt(
            REACTOME_INTERACTIONS.format(organism=ORGANISM["file"].get(
                taxon_identifier, "all_species")),
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

                for primary_interactor_a in primary_accession.get(
                        interactor_a, {interactor_a}):
                    for primary_interactor_b in primary_accession.get(
                            interactor_b, {interactor_b}):
                        if (primary_interactor_a in network
                                and primary_interactor_b in network and
                                primary_interactor_a != primary_interactor_b):
                            network.add_edge(primary_interactor_a,
                                             primary_interactor_b)
                            network.edges[
                                primary_interactor_a,
                                primary_interactor_b]["Reactome"] = 0.5


def add_pathways(network, taxon_identifier=9606):
    for row in download.tabular_txt(REACTOME_PATHWAYS,
                                    delimiter="\t",
                                    usecols=[0, 1, 2]):
        if row[2] == ORGANISM["data"].get(taxon_identifier):
            network.add_node(row[0])
            network.nodes[row[0]]["pathway"] = row[1]


def add_pathway_relations(network):
    for row in download.tabular_txt(REACTOME_PATHWAY_RELATIONS,
                                    delimiter="\t",
                                    usecols=[0, 1]):
        if row[0] in network and row[1] in network:
            network.add_edge(row[0], row[1])
