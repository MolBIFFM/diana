import networkx as nx
import itertools

from pipeline.utilities import download

CORUM_ZIP_ARCHIVE = (
    "https://mips.helmholtz-muenchen.de/corum/download/allComplexes.txt.zip")
CORUM = "allComplexes.txt"


def add_proteins(network, protein_complex_purification_method=[]):
    nodes_to_add = set()
    for row in download.tabular_txt(
            CORUM_ZIP_ARCHIVE,
            file=CORUM,
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
                if (interactor_a in network and interactor_b not in network):
                    nodes_to_add.add(interactor_b)

                elif (interactor_a not in network and interactor_b in network):
                    nodes_to_add.add(interactor_a)
    network.add_nodes_from(nodes_to_add)


def add_interactions(network, protein_complex_purification_method=[]):
    for row in download.tabular_txt(
            CORUM_ZIP_ARCHIVE,
            file=CORUM,
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
                if interactor_a in network and interactor_b in network:
                    network.add_edge(
                        interactor_a,
                        interactor_b,
                    )
                    network.edges[interactor_a, interactor_b, ]["CORUM"] = 0.5
