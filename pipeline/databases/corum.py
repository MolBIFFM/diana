import itertools

from pipeline.utilities import download, uniprot

CORUM_ZIP_ARCHIVE = (
    "https://mips.helmholtz-muenchen.de/corum/download/allComplexes.txt.zip")
CORUM = "allComplexes.txt"


def add_proteins(network, protein_complex_purification_method=[]):
    primary_accessions = uniprot.get_primary_accession()

    nodes_to_add = set()
    for row in download.tabular_txt(
            CORUM_ZIP_ARCHIVE,
            file_from_zip_archive=CORUM,
            delimiter="\t",
            header=0,
            usecols=[
                "subunits(UniProt IDs)",
                "Protein complex purification method",
            ],
    ):
        if not protein_complex_purification_method or any(
                method in protein_complex_purification_method for method in [
                    entry.split("-")[1].lstrip() for entry in
                    row["Protein complex purification method"].split(";")
                ]):
            for interactor_a, interactor_b in itertools.combinations(
                    row["subunits(UniProt IDs)"].split(";"), 2):

                for a in primary_accessions.get(interactor_a, {interactor_a}):
                    for b in primary_accessions.get(interactor_b,
                                                    {interactor_b}):
                        if (a in network and b not in network):
                            nodes_to_add.add(b)

                        elif (a not in network and b in network):
                            nodes_to_add.add(a)

    network.add_nodes_from(nodes_to_add)


def add_interactions(network, protein_complex_purification_method=[]):
    primary_accessions = uniprot.get_primary_accession(network)

    for row in download.tabular_txt(
            CORUM_ZIP_ARCHIVE,
            file_from_zip_archive=CORUM,
            delimiter="\t",
            header=0,
            usecols=[
                "subunits(UniProt IDs)",
                "Protein complex purification method",
            ],
    ):
        if not protein_complex_purification_method or any(
                method in protein_complex_purification_method for method in [
                    entry.split("-")[1].lstrip() for entry in
                    row["Protein complex purification method"].split(";")
                ]):
            for interactor_a, interactor_b in itertools.combinations(
                    row["subunits(UniProt IDs)"].split(";"), 2):

                for a in primary_accessions.get(interactor_a, {interactor_a}):
                    for b in primary_accessions.get(interactor_b,
                                                    {interactor_b}):
                        if a in network and b in network and a != b:
                            network.add_edge(
                                a,
                                b,
                            )
                            network.edges[a, b, ]["CORUM"] = 0.5
