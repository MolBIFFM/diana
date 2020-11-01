from pipeline.download import download


def add_interactions_from_IntAct(self):
    for _, row in download.iterate_tabular_data(
            "ftp://ftp.ebi.ac.uk/pub/databases/intact/current/psimitab/intact.zip",
            delimiter="\t",
            header=0,
            usecols=["#ID(s) interactor A", "ID(s) interactor B"]):
        print(row["#ID(s) interactor A"], row["ID(s) interactor B"])
