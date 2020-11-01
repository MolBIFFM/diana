from pipeline.download import download


def add_interactions_from_STRING(self,
                                 neighborhood=0.0,
                                 neighborhood_transferred=0.0,
                                 fusion=0.0,
                                 cooccurence=0.0,
                                 homology=0.0,
                                 coexpression=0.0,
                                 coexpression_transferred=0.0,
                                 experiments=0.7,
                                 experiments_transferred=0.0,
                                 database=0.0,
                                 database_transferred=0.0,
                                 textmining=0.0,
                                 textmining_transferred=0.0,
                                 combined_score=0.7):
    uniprot_id = {}
    for _, row in download.iterate_tabular_data(
            "ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/by_organism/HUMAN_9606_idmapping.dat.gz",
            delimiter="\t",
            usecols=[0, 1, 2]):
        if row[1] == "STRING" and row[0] in self.nodes:
            uniprot_id[row[2]] = row[0]

    for _, row in download.iterate_tabular_data(
            "https://string-db.org/mapping_files/uniprot/human.uniprot_2_string.2018.tsv.gz",
            usecols=[1, 2]):
        if row[1].split("|")[0] in self.nodes:
            uniprot_id[row[2]] = row[1].split("|")[0]

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
            "combined_score": combined_score
        }.items() if threshold
    }

    for _, row in download.iterate_tabular_data(
            "https://stringdb-static.org/download/protein.links.full.v11.0/9606.protein.links.full.v11.0.txt.gz",
            delimiter=" ",
            header=0,
            usecols=["protein1", "protein2"] + list(thresholds.keys())):
        if (uniprot_id.get(row["protein1"]) and uniprot_id.get(row["protein2"])
                and all(row[column] / 1000 >= thresholds[column]
                        for column in thresholds)):
            self.add_edge(uniprot_id[row["protein1"]],
                          uniprot_id[row["protein2"]])
