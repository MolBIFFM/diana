from pipeline.download import download


def add_interactions_from_BioGRID(self):
    uniprot_id = {}
    for _, row in download.iterate_tabular_data(
            "ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/by_organism/HUMAN_9606_idmapping.dat.gz",
            delimiter="\t",
            usecols=[0, 1, 2]):
        if row[1] == "BioGRID" and row[0] in self.nodes:
            uniprot_id[row[2]] = row[0]
    for _, row in download.iterate_tabular_data(
            "https://downloads.thebiogrid.org/Download/BioGRID/Release-Archive/BIOGRID-4.2.191/BIOGRID-ALL-4.2.191.tab3.zip",
            delimiter="\t",
            header=0,
            usecols=["BioGRID ID Interactor A", "BioGRID ID Interactor B"]):
        print(row["BioGRID ID Interactor A"], row["BioGRID ID Interactor B"])
