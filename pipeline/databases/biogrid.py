import networkx as nx
from pipeline.utilities import download

BIOGRID_ID_MAP_ZIP_ARCHIVE = "https://downloads.thebiogrid.org/Download/BioGRID/Latest-Release/BIOGRID-IDENTIFIERS-LATEST.tab.zip"
BIOGRID_ID_MAP = r"BIOGRID-IDENTIFIERS-[0-9]\.[0-9]\.[0-9]{3}\.tab\.txt"
BIOGRID_ZIP_ARCHIVE = "https://downloads.thebiogrid.org/Download/BioGRID/Latest-Release/BIOGRID-ORGANISM-LATEST.tab3.zip"
BIOGRID = r"BIOGRID-ORGANISM-{organism}-[0-9]\.[0-9]\.[0-9][0-9][0-9]\.tab3\.txt"
BIOGRID_MV_PHYSICAL_ZIP_ARCHIVE = "https://downloads.thebiogrid.org/Download/BioGRID/Latest-Release/BIOGRID-MV-Physical-LATEST.tab3.zip"
BIOGRID_MV_PHYSICAL = r"BIOGRID-MV-Physical-[0-9]\.[0-9]\.[0-9]{3}\.tab3\.txt"

ORGANISM = {
    9606: {
        BIOGRID: "Homo_sapiens"
    },
}


def add_proteins(
    network,
    experimental_system=[],
    experimental_system_type=["physical"],
    taxon_identifier=9606,
    multi_validated_physical=False,
):

    uniprot = {}
    for row in download.tabular_txt(
            BIOGRID_ID_MAP_ZIP_ARCHIVE,
            file=BIOGRID_ID_MAP,
            delimiter="\t",
            header=20,
            usecols=[
                "BIOGRID_ID",
                "IDENTIFIER_VALUE",
                "IDENTIFIER_TYPE",
            ],
    ):
        if row["IDENTIFIER_TYPE"] in ("UNIPROT-ACCESSION", "UNIPROT-ISOFORM"):
            uniprot[int(row["BIOGRID_ID"])] = row["IDENTIFIER_VALUE"]

    for row in download.tabular_txt(
            BIOGRID_MV_PHYSICAL_ZIP_ARCHIVE
            if multi_validated_physical else BIOGRID_ZIP_ARCHIVE,
            file=BIOGRID_MV_PHYSICAL if multi_validated_physical else
            BIOGRID.format(organism=ORGANISM[taxon_identifier][BIOGRID]),
            delimiter="\t",
            header=0,
            usecols=[
                "BioGRID ID Interactor A", "BioGRID ID Interactor B",
                "Experimental System", "Experimental System Type",
                "Organism ID Interactor A", "Organism ID Interactor B",
                "Throughput", "SWISS-PROT Accessions Interactor A",
                "SWISS-PROT Accessions Interactor B"
            ],
    ):

        for interactor_a in row["SWISS-PROT Accessions Interactor A"].split(
                "|"):
            for interactor_b in row[
                    "SWISS-PROT Accessions Interactor B"].split("|"):
                if (interactor_a in network or interactor_b in network
                        and row["Organism ID Interactor A"] !=
                        row["Organism ID Interactor B"] and
                    (not experimental_system
                     or row["Experimental System"] in experimental_system) and
                    (not experimental_system_type
                     or row["Experimental System Type"]
                     in experimental_system_type)):
                    if (interactor_a in network and interactor_b not in network
                            and network.nodes[interactor_a].get("protein")):
                        network.add_node(interactor_b)

                    elif (interactor_a not in network
                          and interactor_b in network
                          and network.nodes[interactor_b].get("protein")):
                        network.add_node(interactor_a)

        if (uniprot.get(row["BioGRID ID Interactor A"]) in network
                or uniprot.get(row["BioGRID ID Interactor B"]) in network
                and uniprot.get(row["BioGRID ID Interactor A"]) != uniprot.get(
                    row["BioGRID ID Interactor B"]) and
            (not experimental_system
             or row["Experimental System"] in experimental_system) and
            (not experimental_system_type
             or row["Experimental System Type"] in experimental_system_type)):
            if (uniprot.get(row["BioGRID ID Interactor A"]) in network
                    and row["BioGRID ID Interactor B"] in uniprot
                    and uniprot[row["BioGRID ID Interactor B"]] not in network
                    and network.nodes[uniprot["BioGRID ID Interactor A"]].get(
                        "protein")):
                network.add_node(uniprot[row["BioGRID ID Interactor B"]])

            elif (row["BioGRID ID Interactor A"] in uniprot
                  and uniprot[row["BioGRID ID Interactor A"]] not in network
                  and uniprot[row["BioGRID ID Interactor B"]] in network
                  and network.nodes[uniprot[
                      row["BioGRID ID Interactor B"]]].get("protein")):
                network.add_node(uniprot[row["BioGRID ID Interactor A"]])


def add_interactions(
    network,
    experimental_system=[],
    experimental_system_type=["physical"],
    taxon_identifier=9606,
    multi_validated_physical=False,
):
    uniprot = {}
    for row in download.tabular_txt(
            BIOGRID_ID_MAP_ZIP_ARCHIVE,
            file=BIOGRID_ID_MAP,
            delimiter="\t",
            header=20,
            usecols=[
                "BIOGRID_ID",
                "IDENTIFIER_VALUE",
                "IDENTIFIER_TYPE",
            ],
    ):
        if row["IDENTIFIER_TYPE"] in ("UNIPROT-ACCESSION", "UNIPROT-ISOFORM"):
            uniprot[int(row["BIOGRID_ID"])] = row["IDENTIFIER_VALUE"]

    for row in download.tabular_txt(
            BIOGRID_MV_PHYSICAL_ZIP_ARCHIVE
            if multi_validated_physical else BIOGRID_ZIP_ARCHIVE,
            file=BIOGRID_MV_PHYSICAL if multi_validated_physical else
            BIOGRID.format(organism=ORGANISM[taxon_identifier][BIOGRID]),
            delimiter="\t",
            header=0,
            usecols=[
                "BioGRID ID Interactor A",
                "BioGRID ID Interactor B",
                "SWISS-PROT Accessions Interactor A",
                "SWISS-PROT Accessions Interactor B",
                "Experimental System",
                "Experimental System Type",
                "Organism ID Interactor A",
                "Organism ID Interactor B",
                "Throughput",
            ],
    ):
        for interactor_a in row["SWISS-PROT Accessions Interactor A"].split(
                "|"):
            for interactor_b in row[
                    "SWISS-PROT Accessions Interactor B"].split("|"):
                if (interactor_a in network and interactor_b in network
                        and interactor_a != interactor_b and
                    (not experimental_system
                     or row["Experimental System"] in experimental_system)
                        and (not experimental_system_type
                             or row["Experimental System Type"]
                             in experimental_system_type)):
                    network.add_edge(
                        interactor_a,
                        interactor_b,
                    )
                    network.edges[interactor_a,
                                  interactor_b, ]["BioGRID"] = 0.5

        if (uniprot.get(row["BioGRID ID Interactor A"]) in network
                and uniprot.get(row["BioGRID ID Interactor B"]) in network
                and uniprot[row["BioGRID ID Interactor A"]] !=
                uniprot[row["BioGRID ID Interactor B"]]
                and (not experimental_system
                     or row["Experimental System"] in experimental_system) and
            (not experimental_system_type
             or row["Experimental System Type"] in experimental_system_type)):
            network.add_edge(
                uniprot[row["BioGRID ID Interactor A"]],
                uniprot[row["BioGRID ID Interactor B"]],
            )
            network.edges[
                uniprot[row["BioGRID ID Interactor A"]],
                uniprot[row["BioGRID ID Interactor B"]], ]["BioGRID"] = 0.5
