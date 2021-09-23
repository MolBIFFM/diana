from pipeline.utilities import download, uniprot

BIOGRID_ID_MAP_ZIP_ARCHIVE = "https://downloads.thebiogrid.org/Download/BioGRID/Latest-Release/BIOGRID-IDENTIFIERS-LATEST.tab.zip"
BIOGRID_ID_MAP = r"BIOGRID-IDENTIFIERS-[0-9]\.[0-9]\.[0-9]{3}\.tab\.txt"
BIOGRID_ZIP_ARCHIVE = "https://downloads.thebiogrid.org/Download/BioGRID/Latest-Release/BIOGRID-ORGANISM-LATEST.tab3.zip"
BIOGRID = r"BIOGRID-ORGANISM-{organism}-[0-9]\.[0-9]\.[0-9][0-9][0-9]\.tab3\.txt"
BIOGRID_MV_PHYSICAL_ZIP_ARCHIVE = "https://downloads.thebiogrid.org/Download/BioGRID/Latest-Release/BIOGRID-MV-Physical-LATEST.tab3.zip"
BIOGRID_MV_PHYSICAL = r"BIOGRID-MV-Physical-[0-9]\.[0-9]\.[0-9]{3}\.tab3\.txt"

ORGANISM = {
    3055: "Chlamydomonas_reinhardtii",
    3702: "Arabidopsis_thaliana_Columbia",
    3847: "Glycine_max",
    3988: "Ricinus_communis",
    4081: "Solanum_lycopersicum",
    4098: "Nicotiana_tomentosiformis",
    4113: "Solanum_tuberosum",
    4558: "Sorghum_bicolor",
    4577: "Zea_mays",
    6239: "Caenorhabditis_elegans",
    7227: "Drosophila_melanogaster",
    7460: "Apis_mellifera",
    7668: "Strongylocentrotus_purpuratus",
    7955: "Danio_rerio",
    8355: "Xenopus_laevis",
    8839: "Anas_platyrhynchos",
    9031: "Gallus_gallus",
    9103: "Meleagris_gallopavo",
    9544: "Macaca_mulatta",
    9598: "Pan_troglodytes",
    9606: "Homo_sapiens",
    9615: "Canis_familiaris",
    9685: "Felis_catus",
    9796: "Equus_caballus",
    9823: "Sus_scrofa",
    9940: "Ovis_aries",
    9986: "Oryctolagus_cuniculus",
    10029: "Cricetulus_griseus",
    10090: "Mus_musculus",
    10116: "Rattus_norvegicus",
    10141: "Cavia_porcellus",
    10245: "Vaccinia_Virus",
    10298: "Human_Herpesvirus_1",
    10310: "Human_Herpesvirus_2",
    10335: "Human_Herpesvirus_3",
    10359: "Human_Herpesvirus_5",
    10372: "Human_Herpesvirus_7",
    10376: "Human_Herpesvirus_4",
    10600: "Human_papillomavirus_6b",
    10620: "Human_papillomavirus_7",
    10621: "Human_papillomavirus_9",
    10633: "Simian_Virus_40",
    11103: "Hepatitus_C_Virus",
    11676: "Human_Immunodeficiency_Virus_1",
    11709: "Human_Immunodeficiency_Virus_2",
    11723: "Simian_Immunodeficiency_Virus",
    12242: "Tobacco_Mosaic_Virus",
    13616: "Monodelphis_domestica",
    29760: "Vitis_vinifera",
    32603: "Human_Herpesvirus_6A",
    32604: "Human_Herpesvirus_6B",
    36329: "Plasmodium_falciparum_3D7",
    37296: "Human_Herpesvirus_8",
    39947: "Oryza_sativa_Japonica",
    44689: "Dictyostelium_discoideum",
    60711: "Chlorocebus_sabaeus",
    83332: "Mycobacterium_tuberculosis_H37Rv",
    83333: "Escherichia_coli_K12",
    88036: "Selaginella_moellendorffii",
    121224: "Pediculus_humanus",
    171101: "Streptococcus_pneumoniae_ATCCBAA255",
    180454: "Anopheles_gambiae_PEST",
    224308: "Bacillus_subtilis_168",
    227321: "Emericella_nidulans_FGSC_A4",
    237561: "Candida_albicans_SC5314",
    237631: "Ustilago_maydis_521",
    284812: "Schizosaccharomyces_pombe_972h",
    316407: "Escherichia_coli_K12_W3110",
    333759: "Human_papillomavirus_10",
    333760: "Human_papillomavirus_16",
    333763: "Human_papillomavirus_32",
    333923: "Human_papillomavirus_5",
    347515: "Leishmania_major_Friedlin",
    367110: "Neurospora_crassa_OR74A",
    511145: "Escherichia_coli_K12_MG1655",
    559292: "Saccharomyces_cerevisiae_S288c",
    595496: "Escherichia_coli_K12_MC4100_BW2952",
    694009: "Severe_acute_respiratory_syndrome_coronavirus",
    1335626: "Middle-East_Respiratory_Syndrome-related_Coronavirus",
    2697049: "Severe_acute_respiratory_syndrome_coronavirus_2",
}


def add_proteins(
    network,
    experimental_system=[],
    experimental_system_type=[],
    multi_validated_physical=False,
    taxon_identifier=9606,
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

    primary_accession = uniprot.get_primary_accession(taxon_identifier)

    nodes_to_add = set()
    for row in download.tabular_txt(
            BIOGRID_MV_PHYSICAL_ZIP_ARCHIVE
            if multi_validated_physical else BIOGRID_ZIP_ARCHIVE,
            file_from_zip_archive=BIOGRID_MV_PHYSICAL
            if multi_validated_physical else BIOGRID.format(
                organism=ORGANISM[taxon_identifier]),
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
                for int_a in primary_accession.get(interactor_a,
                                                   {interactor_a}):
                    for int_b in primary_accession.get(interactor_b,
                                                       {interactor_b}):
                        if ((not experimental_system or
                             row["Experimental System"] in experimental_system)
                                and (not experimental_system_type
                                     or row["Experimental System Type"]
                                     in experimental_system_type)):
                            if (int_a in network and int_b not in network):
                                nodes_to_add.add(int_b)

                            elif (int_a not in network and int_b in network):
                                nodes_to_add.add(int_a)

        for interactor_a in uniprot_id_map.get(row["BioGRID ID Interactor A"],
                                               {}):
            for interactor_b in uniprot_id_map.get(
                    row["BioGRID ID Interactor B"], {}):
                for int_a in primary_accession.get(interactor_a,
                                                   {interactor_a}):
                    for int_b in primary_accession.get(interactor_b,
                                                       {interactor_b}):
                        if (int_a != int_b and int_a in network
                                or int_b in network and
                            (not experimental_system or
                             row["Experimental System"] in experimental_system)
                                and (not experimental_system_type
                                     or row["Experimental System Type"]
                                     in experimental_system_type)):
                            if (int_a in network and int_b not in network):
                                nodes_to_add.add(int_b)

                            elif (int_a not in network and int_b in network):
                                nodes_to_add.add(int_b)

    network.add_nodes_from(nodes_to_add)


def add_interactions(
    network,
    experimental_system=[],
    experimental_system_type=[],
    multi_validated_physical=False,
    taxon_identifier=9606,
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

    primary_accession = uniprot.get_primary_accession(taxon_identifier,
                                                      network)

    for row in download.tabular_txt(
            BIOGRID_MV_PHYSICAL_ZIP_ARCHIVE
            if multi_validated_physical else BIOGRID_ZIP_ARCHIVE,
            file_from_zip_archive=BIOGRID_MV_PHYSICAL
            if multi_validated_physical else BIOGRID.format(
                organism=ORGANISM[taxon_identifier]),
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
                for int_a in primary_accession.get(interactor_a,
                                                   {interactor_a}):
                    for int_b in primary_accession.get(interactor_b,
                                                       {interactor_b}):
                        if (int_a in network and int_b in network
                                and int_a != int_b and
                            (not experimental_system or
                             row["Experimental System"] in experimental_system)
                                and (not experimental_system_type
                                     or row["Experimental System Type"]
                                     in experimental_system_type)):
                            network.add_edge(int_a, int_b)
                            network.edges[int_a, int_b]["BioGRID"] = 0.5

        for interactor_a in uniprot_id_map.get(row["BioGRID ID Interactor A"],
                                               {}):
            for interactor_b in uniprot_id_map.get(
                    row["BioGRID ID Interactor B"], {}):
                for int_a in primary_accession.get(interactor_a,
                                                   {interactor_a}):
                    for int_b in primary_accession.get(interactor_b,
                                                       {interactor_b}):
                        if ((not experimental_system or
                             row["Experimental System"] in experimental_system)
                                and (not experimental_system_type
                                     or row["Experimental System Type"]
                                     in experimental_system_type)
                                and int_a != int_b):
                            network.add_edge(int_a, int_b)
                            network.edges[int_a, int_b]["BioGRID"] = 0.5
