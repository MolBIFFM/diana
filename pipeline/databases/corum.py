import itertools

from pipeline.utilities import download, uniprot

CORUM_ZIP_ARCHIVE = (
    "https://mips.helmholtz-muenchen.de/corum/download/allComplexes.txt.zip")
CORUM = "allComplexes.txt"


def add_proteins(network,
                 protein_complex_purification_method=[],
                 taxon_identifier=9606):
    primary_accessions = uniprot.get_primary_accession(taxon_identifier)

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

                for int_a in primary_accessions.get(interactor_a,
                                                    {interactor_a}):
                    for int_b in primary_accessions.get(
                            interactor_b, {interactor_b}):
                        if (int_a in network and int_b not in network):
                            nodes_to_add.add(int_b)

                        elif (int_a not in network and int_b in network):
                            nodes_to_add.add(int_a)

    network.add_nodes_from(nodes_to_add)


def add_interactions(network,
                     protein_complex_purification_method=[],
                     taxon_identifier=9606):
    primary_accessions = uniprot.get_primary_accession(taxon_identifier,
                                                       network)

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

                for int_a in primary_accessions.get(interactor_a,
                                                    {interactor_a}):
                    for int_b in primary_accessions.get(
                            interactor_b, {interactor_b}):
                        if int_a in network and int_b in network and int_a != int_b:
                            network.add_edge(
                                int_a,
                                int_b,
                            )
                            network.edges[int_a, int_b, ]["CORUM"] = 0.5
