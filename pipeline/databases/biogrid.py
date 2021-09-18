from pipeline.utilities import download, uniprot

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
    uniprot_id_map = uniprot.get_id_map("BioGRID", taxon_identifier, set(),
                                        int)

    for row in download.tabular_txt(
            BIOGRID_ID_MAP_ZIP_ARCHIVE,
            file_from_zip_archive=BIOGRID_ID_MAP,
            delimiter="\t",
            header=20,
            usecols=[
                "BIOGRID_ID",
                "IDENTIFIER_VALUE",
                "IDENTIFIER_TYPE",
            ],
    ):
        if row["IDENTIFIER_TYPE"] in ("UNIPROT-ACCESSION", "UNIPROT-ISOFORM"):
            if int(row["BIOGRID_ID"]) not in uniprot_id_map:
                uniprot_id_map[int(row["BIOGRID_ID"])] = set()
            uniprot_id_map[int(row["BIOGRID_ID"])].add(row["IDENTIFIER_VALUE"])

    primary_accession = uniprot.get_primary_accession()

    nodes_to_add = set()
    for row in download.tabular_txt(
            BIOGRID_MV_PHYSICAL_ZIP_ARCHIVE
            if multi_validated_physical else BIOGRID_ZIP_ARCHIVE,
            file_from_zip_archive=BIOGRID_MV_PHYSICAL
            if multi_validated_physical else BIOGRID.format(
                organism=ORGANISM[taxon_identifier][BIOGRID]),
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
                for a in primary_accession.get(interactor_a, {interactor_a}):
                    for b in primary_accession.get(interactor_b,
                                                   {interactor_b}):
                        if (a in network or b in network and a != b and
                            (not experimental_system or
                             row["Experimental System"] in experimental_system)
                                and (not experimental_system_type
                                     or row["Experimental System Type"]
                                     in experimental_system_type)):
                            if (a in network and b not in network):
                                nodes_to_add.add(b)

                            elif (a not in network and b in network):
                                nodes_to_add.add(a)

        for protein_a in uniprot_id_map.get(row["BioGRID ID Interactor A"],
                                            {}):
            for protein_b in uniprot_id_map.get(row["BioGRID ID Interactor B"],
                                                {}):
                for a in primary_accession.get(protein_a, {protein_a}):
                    for b in primary_accession.get(protein_b, {protein_b}):
                        if (a in network or b in network and a != b and
                            (not experimental_system or
                             row["Experimental System"] in experimental_system)
                                and (not experimental_system_type
                                     or row["Experimental System Type"]
                                     in experimental_system_type)):
                            if (a in network and b and b not in network):
                                nodes_to_add.add(b)

                            elif (a and a not in network and b in network):
                                nodes_to_add.add(a)

    network.add_nodes_from(nodes_to_add)


def add_interactions(
    network,
    experimental_system=[],
    experimental_system_type=["physical"],
    taxon_identifier=9606,
    multi_validated_physical=False,
):
    uniprot_id_map = uniprot.get_id_map("BioGRID", taxon_identifier, network,
                                        int)

    test = {}
    for row in download.tabular_txt(
            BIOGRID_ID_MAP_ZIP_ARCHIVE,
            file_from_zip_archive=BIOGRID_ID_MAP,
            delimiter="\t",
            header=20,
            usecols=[
                "BIOGRID_ID",
                "IDENTIFIER_VALUE",
                "IDENTIFIER_TYPE",
            ],
    ):

        if row["IDENTIFIER_TYPE"] in ("UNIPROT-ACCESSION", "UNIPROT-ISOFORM"
                                      ) and row["IDENTIFIER_VALUE"] in network:
            if int(row["BIOGRID_ID"]) not in uniprot_id_map:
                uniprot_id_map[int(row["BIOGRID_ID"])] = set()
            uniprot_id_map[int(row["BIOGRID_ID"])].add(row["IDENTIFIER_VALUE"])

    primary_accession = uniprot.get_primary_accession(network)

    for row in download.tabular_txt(
            BIOGRID_MV_PHYSICAL_ZIP_ARCHIVE
            if multi_validated_physical else BIOGRID_ZIP_ARCHIVE,
            file_from_zip_archive=BIOGRID_MV_PHYSICAL
            if multi_validated_physical else BIOGRID.format(
                organism=ORGANISM[taxon_identifier][BIOGRID]),
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
                for a in primary_accession.get(interactor_a, {interactor_a}):
                    for b in primary_accession.get(interactor_b,
                                                   {interactor_b}):
                        if (a in network and b in network and a != b and
                            (not experimental_system or
                             row["Experimental System"] in experimental_system)
                                and (not experimental_system_type
                                     or row["Experimental System Type"]
                                     in experimental_system_type)):
                            network.add_edge(
                                a,
                                b,
                            )
                            network.edges[a, b, ]["BioGRID"] = 0.5

        for protein_a in uniprot_id_map.get(row["BioGRID ID Interactor A"],
                                            {}):
            for protein_b in uniprot_id_map.get(row["BioGRID ID Interactor B"],
                                                {}):
                for a in primary_accession.get(protein_a, {protein_a}):
                    for b in primary_accession.get(protein_b, {protein_b}):
                        if (a != b and
                            (not experimental_system or
                             row["Experimental System"] in experimental_system)
                                and (not experimental_system_type
                                     or row["Experimental System Type"]
                                     in experimental_system_type)):
                            network.add_edge(
                                a,
                                b,
                            )
                            network.edges[a, b, ]["BioGRID"] = 0.5
