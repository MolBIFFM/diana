from pipeline.utilities import download, uniprot

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
    uniprot_id_map = uniprot.get_id_map("STRING", taxon_identifier)

    for row in download.tabular_txt(
            STRING_ID_MAP.format(taxon_identifier=taxon_identifier),
            delimiter="\t",
            skiprows=1,
            usecols=[0, 1, 2],
    ):
        if "BLAST_UniProt_AC" in row[2].split():
            if row[0] not in uniprot_id_map:
                uniprot_id_map[row[0]] = set()
            uniprot_id_map[row[0]].add(row[1])

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

    primary_accession = uniprot.get_primary_accession()

    nodes_to_add = set()
    for row in download.tabular_txt(
            STRING_PHYSICAL.format(taxon_identifier=taxon_identifier)
            if physical else STRING.format(taxon_identifier=taxon_identifier),
            delimiter=" ",
            header=0,
            usecols=["protein1", "protein2"] + list(thresholds.keys()),
    ):
        for protein_a in uniprot_id_map.get(row["protein1"], {}):
            for protein_b in uniprot_id_map.get(row["protein2"], {}):
                for a in primary_accession.get(protein_a, {protein_a}):
                    for b in primary_accession.get(protein_b, {protein_b}):
                        if (a in network or b in network and a != b and all(
                                row[column] / 1000 >= thresholds[column]
                                for column in thresholds)):
                            if (a in network and b and b not in network):
                                nodes_to_add.add(b)

                            elif (a and a not in network and b in network):
                                nodes_to_add.add(a)

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
    uniprot_id_map = uniprot.get_id_map("STRING", taxon_identifier, network)

    for row in download.tabular_txt(
            STRING_ID_MAP.format(taxon_identifier=taxon_identifier),
            delimiter="\t",
            skiprows=1,
            usecols=[0, 1, 2],
    ):
        if "BLAST_UniProt_AC" in row[2].split() and row[1] in network:
            if row[0] not in uniprot_id_map:
                uniprot_id_map[row[0]] = set()
            uniprot_id_map[row[0]].add(row[1])

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

    primary_accession = uniprot.get_primary_accession(network)

    for row in download.tabular_txt(
            STRING_PHYSICAL.format(taxon_identifier=taxon_identifier)
            if physical else STRING.format(taxon_identifier=taxon_identifier),
            delimiter=" ",
            header=0,
            usecols=["protein1", "protein2"] + list(thresholds.keys()),
    ):
        for protein_a in uniprot_id_map.get(row["protein1"], {}):
            for protein_b in uniprot_id_map.get(row["protein2"], {}):
                for a in primary_accession.get(protein_a, {protein_a}):
                    for b in primary_accession.get(protein_b, {protein_b}):
                        if (a != b and all(
                                row[column] / 1000 >= thresholds[column]
                                for column in thresholds)):

                            if network.has_edge(a, b):
                                network.edges[a, b]["STRING"] = max(
                                    row["combined_score"] / 1000,
                                    network.edges[a, b].get("STRING", 0.0),
                                )
                            else:
                                network.add_edge(a, b)
                                network.edges[a, b]["STRING"] = (
                                    row["combined_score"] / 1000)
