import networkx as nx

from pipeline.utilities import download

STRING_ID_MAP = "https://stringdb-static.org/download/protein.aliases.v11.0/{taxon_identifier}.protein.aliases.v11.0.txt.gz"
STRING = "https://stringdb-static.org/download/protein.links.full.v11.0/{taxon_identifier}.protein.links.full.v11.0.txt.gz"
STRING_PHYSICAL = "https://stringdb-static.org/download/protein.physical.links.full.v11.0/{taxon_identifier}.protein.links.full.v11.0.txt.gz"


def add_proteins(
    network,
    neighborhood=0.0,
    neighborhood_transferred=0.0,
    fusion=0.0,
    cooccurence=0.0,
    homology=0.0,
    coexpression=0.0,
    coexpression_transferred=0.0,
    experiments=0.0,
    experiments_transferred=0.0,
    database=0.0,
    database_transferred=0.0,
    textmining=0.0,
    textmining_transferred=0.0,
    combined_score=0.0,
    taxon_identifier=9606,
    physical=False,
):
    uniprot = {}
    for row in download.tabular_txt(
            STRING_ID_MAP.format(taxon_identifier=taxon_identifier),
            delimiter="\t",
            skiprows=1,
            usecols=[0, 1, 2],
    ):
        if "BLAST_UniProt_AC" in row[2].split():
            uniprot[row[0]] = row[1]

    thresholds = {
        column: threshold
        for column, threshold in {
            "neighborhood": neighborhood,
            "neighborhood_transferred": neighborhood_transferred,
            "fusion": fusion,
            "cooccurence": cooccurence,
            "homology": homology,
            "coexpression": coexpression,
            "coexpression_transferred": coexpression_transferred,
            "experiments": experiments,
            "experiments_transferred": experiments_transferred,
            "database": database,
            "database_transferred": database_transferred,
            "textmining": textmining,
            "textmining_transferred": textmining_transferred,
        }.items() if threshold
    }
    thresholds["combined_score"] = combined_score

    nodes_to_add = set()
    for row in download.tabular_txt(
            STRING_PHYSICAL.format(taxon_identifier=taxon_identifier)
            if physical else STRING.format(taxon_identifier=taxon_identifier),
            delimiter=" ",
            header=0,
            usecols=["protein1", "protein2"] + list(thresholds.keys()),
    ):
        if (row["protein1"] in uniprot and row["protein2"] in uniprot
                and uniprot[row["protein1"]] != uniprot[row["protein2"]]
                and all(row[column] / 1000 >= thresholds[column]
                        for column in thresholds)):
            if (uniprot[row["protein1"]] in network
                    and uniprot[row["protein2"]] not in network):
                nodes_to_add.add(uniprot[row["protein2"]])

            elif (uniprot[row["protein1"]] not in network
                  and uniprot[row["protein2"]] in network):
                nodes_to_add.add(uniprot[row["protein1"]])
    network.add_nodes_from(nodes_to_add)


def add_interactions(
    network,
    neighborhood=0.0,
    neighborhood_transferred=0.0,
    fusion=0.0,
    cooccurence=0.0,
    homology=0.0,
    coexpression=0.0,
    coexpression_transferred=0.0,
    experiments=0.0,
    experiments_transferred=0.0,
    database=0.0,
    database_transferred=0.0,
    textmining=0.0,
    textmining_transferred=0.0,
    combined_score=0.0,
    taxon_identifier=9606,
    physical=False,
):

    uniprot = {}
    for row in download.tabular_txt(
            STRING_ID_MAP.format(taxon_identifier=taxon_identifier),
            delimiter="\t",
            skiprows=1,
            usecols=[0, 1, 2],
    ):
        if "BLAST_UniProt_AC" in row[2].split() and row[1] in network:
            uniprot[row[0]] = row[1]

    thresholds = {
        column: threshold
        for column, threshold in {
            "neighborhood": neighborhood,
            "neighborhood_transferred": neighborhood_transferred,
            "fusion": fusion,
            "cooccurence": cooccurence,
            "homology": homology,
            "coexpression": coexpression,
            "coexpression_transferred": coexpression_transferred,
            "experiments": experiments,
            "experiments_transferred": experiments_transferred,
            "database": database,
            "database_transferred": database_transferred,
            "textmining": textmining,
            "textmining_transferred": textmining_transferred,
        }.items() if threshold
    }
    thresholds["combined_score"] = combined_score

    for row in download.tabular_txt(
            STRING_PHYSICAL.format(taxon_identifier=taxon_identifier)
            if physical else STRING.format(taxon_identifier=taxon_identifier),
            delimiter=" ",
            header=0,
            usecols=["protein1", "protein2"] + list(thresholds.keys()),
    ):
        if (row["protein1"] in uniprot and row["protein2"] in uniprot
                and uniprot[row["protein1"]] != uniprot[row["protein2"]]
                and all(row[column] / 1000 >= thresholds[column]
                        for column in thresholds)):
            if network.has_edge(uniprot[row["protein1"]],
                                uniprot[row["protein2"]]):
                network.edges[uniprot[row["protein1"]],
                              uniprot[row["protein2"]]]["STRING"] = max(
                                  row["combined_score"] / 1000,
                                  network.edges[uniprot[row["protein1"]],
                                                uniprot[row["protein2"]]].get(
                                                    "STRING", 0.0),
                              )
            else:
                network.add_edge(uniprot[row["protein1"]],
                                 uniprot[row["protein2"]])
                network.edges[uniprot[row["protein1"]],
                              uniprot[row["protein2"]]]["STRING"] = (
                                  row["combined_score"] / 1000)
