import itertools

from uniprot import uniprot
from download import download

CORUM_ZIP_ARCHIVE = "https://mips.helmholtz-muenchen.de/corum/download/allComplexes.txt.zip"
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

                for primary_interactor_a in primary_accessions.get(
                        interactor_a, set()):
                    for primary_interactor_b in primary_accessions.get(
                            interactor_b, set()):
                        if (primary_interactor_a in network
                                and primary_interactor_b not in network):
                            nodes_to_add.add(primary_interactor_b)

                        elif (primary_interactor_a not in network
                              and primary_interactor_b in network):
                            nodes_to_add.add(primary_interactor_a)

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

                for primary_interactor_a in primary_accessions.get(
                        interactor_a, set()):
                    for primary_interactor_b in primary_accessions.get(
                            interactor_b, set()):
                        if (primary_interactor_a in network
                                and primary_interactor_b in network and
                                primary_interactor_a != primary_interactor_b):
                            network.add_edge(
                                primary_interactor_a,
                                primary_interactor_b,
                            )
                            network.edges[
                                primary_interactor_a,
                                primary_interactor_b, ]["CORUM"] = 1.0
