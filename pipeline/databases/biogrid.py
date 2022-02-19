from databases import uniprot
from download import download

BIOGRID_ID_MAP_ZIP_ARCHIVE = "https://downloads.thebiogrid.org/Download/BioGRID/Latest-Release/BIOGRID-IDENTIFIERS-LATEST.tab.zip"
BIOGRID_ID_MAP = r"BIOGRID-IDENTIFIERS-[0-9]\.[0-9]\.[0-9]{3}\.tab\.txt"
BIOGRID_ZIP_ARCHIVE = "https://downloads.thebiogrid.org/Download/BioGRID/Latest-Release/BIOGRID-ORGANISM-LATEST.tab3.zip"
BIOGRID_INTERACTIONS = r"BIOGRID-ORGANISM-{organism}-[0-9]\.[0-9]\.[0-9][0-9][0-9]\.tab3\.txt"
BIOGRID_MV_PHYSICAL_ZIP_ARCHIVE = "https://downloads.thebiogrid.org/Download/BioGRID/Latest-Release/BIOGRID-MV-Physical-LATEST.tab3.zip"
BIOGRID_MV_PHYSICAL = r"BIOGRID-MV-Physical-[0-9]\.[0-9]\.[0-9]{3}\.tab3\.txt"

ORGANISM = {"file": {9606: "Homo_sapiens"}}


def add_proteins(
    network,
    experimental_system=[],
    experimental_system_type=[],
    interaction_throughput=[],
    multi_validated_physical=False,
    taxon_identifier=9606,
):
    primary_accession = uniprot.get_primary_accession(taxon_identifier)

    nodes_to_add = set()
    for row in download.tabular_txt(
            BIOGRID_MV_PHYSICAL_ZIP_ARCHIVE
            if multi_validated_physical else BIOGRID_ZIP_ARCHIVE,
            file_from_zip_archive=BIOGRID_MV_PHYSICAL
            if multi_validated_physical else BIOGRID_INTERACTIONS.format(
                organism=ORGANISM["file"][taxon_identifier]),
            delimiter="\t",
            header=0,
            usecols=[
                "Experimental System", "Experimental System Type",
                "Organism ID Interactor A", "Organism ID Interactor B",
                "Throughput"
                "SWISS-PROT Accessions Interactor A",
                "SWISS-PROT Accessions Interactor B"
            ],
    ):
        if (row["SWISS-PROT Accessions Interactor A"] == "-"
                or row["SWISS-PROT Accessions Interactor B"] == "-"):
            continue

        if ((not experimental_system
             or row["Experimental System"] in experimental_system) and
            (not experimental_system_type
             or row["Experimental System Type"] in experimental_system_type)
                and row["Organism ID Interactor A"] == taxon_identifier
                and row["Organism ID Interactor B"] == taxon_identifier
                and (not interaction_throughput
                     or any(it in interaction_throughput
                            for it in row["Throughput"].split("|")))):
            for interactor_a in row[
                    "SWISS-PROT Accessions Interactor A"].split("|"):

                if "-" in interactor_a and not interactor_a.split(
                        "-")[1].isnumeric():
                    interactor_a = interactor_a.split("-")[0]

                for interactor_b in row[
                        "SWISS-PROT Accessions Interactor B"].split("|"):

                    if "-" in interactor_b and not interactor_b.split(
                            "-")[1].isnumeric():
                        interactor_b = interactor_b.split("-")[0]

                    for primary_interactor_a in primary_accession.get(
                            interactor_a, {interactor_a}):
                        for primary_interactor_b in primary_accession.get(
                                interactor_b, {interactor_b}):
                            if (primary_interactor_a in network
                                    and primary_interactor_b not in network):
                                nodes_to_add.add(primary_interactor_b)

                            elif (primary_interactor_a not in network
                                  and primary_interactor_b in network):
                                nodes_to_add.add(primary_interactor_a)

    network.add_nodes_from(nodes_to_add)


def add_protein_protein_interactions(
    network,
    experimental_system=[],
    experimental_system_type=[],
    interaction_throughput=[],
    multi_validated_physical=False,
    taxon_identifier=9606,
):
    primary_accession = uniprot.get_primary_accession(taxon_identifier,
                                                      network)

    for row in download.tabular_txt(
            BIOGRID_MV_PHYSICAL_ZIP_ARCHIVE
            if multi_validated_physical else BIOGRID_ZIP_ARCHIVE,
            file_from_zip_archive=BIOGRID_MV_PHYSICAL
            if multi_validated_physical else BIOGRID_INTERACTIONS.format(
                organism=ORGANISM["file"][taxon_identifier]),
            delimiter="\t",
            header=0,
            usecols=[
                "Experimental System", "Experimental System Type",
                "Throughput", "SWISS-PROT Accessions Interactor A",
                "SWISS-PROT Accessions Interactor B"
            ],
    ):
        if (row["SWISS-PROT Accessions Interactor A"] == "-"
                or row["SWISS-PROT Accessions Interactor B"] == "-"):
            continue

        if ((not experimental_system
             or row["Experimental System"] in experimental_system) and
            (not experimental_system_type
             or row["Experimental System Type"] in experimental_system_type)
                and (not interaction_throughput
                     or any(it in interaction_throughput
                            for it in row["Throughput"].split("|")))):
            for interactor_a in row[
                    "SWISS-PROT Accessions Interactor A"].split("|"):

                if "-" in interactor_a and not interactor_a.split(
                        "-")[1].isnumeric():
                    interactor_a = interactor_a.split("-")[0]

                for interactor_b in row[
                        "SWISS-PROT Accessions Interactor B"].split("|"):

                    if "-" in interactor_b and not interactor_b.split(
                            "-")[1].isnumeric():
                        interactor_b = interactor_b.split("-")[0]

                    for primary_interactor_a in primary_accession.get(
                            interactor_a, {interactor_a}):
                        for primary_interactor_b in primary_accession.get(
                                interactor_b, {interactor_b}):
                            if (primary_interactor_a in network
                                    and primary_interactor_b in network
                                    and primary_interactor_a !=
                                    primary_interactor_b):
                                network.add_edge(primary_interactor_a,
                                                 primary_interactor_b)
                                network.edges[
                                    primary_interactor_a,
                                    primary_interactor_b]["BioGRID"] = 1.0
