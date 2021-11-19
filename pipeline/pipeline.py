import argparse
import concurrent.futures
import itertools
import json
import logging
import os

import networkx as nx

from cytoscape import styles
from databases import biogrid, corum, intact, mint, string
from interface import algorithm, combination, conversion, extraction
from networks import protein_interaction_network


def process(configuration_file, log=False):
    if log:
        logger = logging.getLogger("root")
        logging.basicConfig(filename="{}.log".format(
            os.path.splitext(os.path.basename(configuration_file))[0]),
                            filemode="w",
                            level=logging.DEBUG,
                            format="%(message)s")

    with open(configuration_file) as configuration:
        configurations = json.load(configuration)

    for i, configuration in enumerate(configurations, start=1):
        network = nx.Graph()

        for entry in configuration.get("genes", {}):
            if "file" in entry and "accession column" in entry:
                protein_interaction_network.add_genes_from_table(
                    network,
                    file_name=entry["file"],
                    gene_accession_column=entry["accession column"],
                    gene_accession_format=extraction.EXTRACT.get(
                        entry.get("accession format"), lambda entry: [entry]),
                    sheet_name=entry.get("sheet", 0),
                    header=entry.get("header", 1) - 1,
                    taxon_identifier=entry.get("taxon identifier", 9606),
                )

            elif "accessions" in entry:
                protein_interaction_network.add_genes_from(
                    network,
                    genes=entry["accessions"],
                    taxon_identifier=entry.get("taxon identifier", 9606),
                )

        for entry in configuration.get("proteins", {}):
            if "file" in entry and "accession column" in entry:
                protein_interaction_network.add_proteins_from_table(
                    network,
                    file_name=entry["file"],
                    protein_accession_column=entry["accession column"],
                    protein_accession_format=extraction.EXTRACT.get(
                        entry.get("accession format"), lambda entry: [entry]),
                    time=entry.get("time", 0),
                    modification=entry.get("post-translational modification",
                                           ""),
                    position_column=entry.get("position column", ""),
                    position_format=extraction.EXTRACT.get(
                        entry.get("position format"), lambda entry: [entry]),
                    replicates=entry.get("replicate columns", []),
                    sheet_name=entry.get("sheet", 0),
                    header=entry.get("header", 1) - 1,
                    num_sites=entry.get("sites", 1000),
                    num_replicates=entry.get("replicates", 1),
                    replicate_combination=combination.REPLICATE_COMBINATION.
                    get(
                        entry.get("combine replicates", "mean"),
                        combination.REPLICATE_COMBINATION["mean"],
                    ),
                    measurement_conversion=conversion.LOG_BASE[entry.get(
                        "log base")],
                )

            elif "accessions" in entry:
                protein_interaction_network.add_proteins_from(
                    network, proteins=entry["accessions"])

        for entry in configuration.get("networks", []):
            network = nx.algorithms.operators.binary.compose(
                network, nx.readwrite.graphml.read_graphml(entry))

        if "protein-protein interactions" in configuration:
            k = 0
            while any(configuration["protein-protein interactions"].get(
                    database, {}).get("neighbors", 0) > k for database in
                      {"BioGRID", "CORUM", "IntAct", "Reactome", "STRING"}):
                if "BioGRID" in configuration[
                        "protein-protein interactions"] and configuration[
                            "protein-protein interactions"]["BioGRID"].get(
                                "neighbors", 0) > k:
                    biogrid.add_proteins(
                        network,
                        interaction_throughput=configuration[
                            "protein-protein interactions"]["BioGRID"].get(
                                "interaction throughput", []),
                        experimental_system=configuration[
                            "protein-protein interactions"]["BioGRID"].get(
                                "experimental system", []),
                        experimental_system_type=configuration[
                            "protein-protein interactions"]["BioGRID"].get(
                                "experimental system type", []),
                        multi_validated_physical=configuration[
                            "protein-protein interactions"]["BioGRID"].get(
                                "multi-validated physical", False),
                        taxon_identifier=configuration[
                            "protein-protein interactions"]["BioGRID"].get(
                                "taxon identifier", 9606),
                    )

                if "CORUM" in configuration[
                        "protein-protein interactions"] and configuration[
                            "protein-protein interactions"]["CORUM"].get(
                                "neighbors", 0) > k:
                    corum.add_proteins(
                        network,
                        protein_complex_purification_method=configuration[
                            "protein-protein interactions"]["CORUM"].get(
                                "protein complex purification method",
                                [],
                            ),
                        taxon_identifier=configuration[
                            "protein-protein interactions"]["CORUM"].get(
                                "taxon identifier", 9606),
                    )

                if "IntAct" in configuration[
                        "protein-protein interactions"] and configuration[
                            "protein-protein interactions"]["IntAct"].get(
                                "neighbors", 0) > k:
                    intact.add_proteins(
                        network,
                        interaction_detection_methods=configuration[
                            "protein-protein interactions"]["IntAct"].get(
                                "interaction detection methods", []),
                        interaction_types=configuration[
                            "protein-protein interactions"]["IntAct"].get(
                                "interaction types", []),
                        mi_score=configuration["protein-protein interactions"]
                        ["IntAct"].get("MI score", 0.0),
                        taxon_identifier=configuration[
                            "protein-protein interactions"]["IntAct"].get(
                                "taxon identifier", 9606),
                    )

                if "MINT" in configuration[
                        "protein-protein interactions"] and configuration[
                            "protein-protein interactions"]["MINT"].get(
                                "neighbors", 0) > k:
                    mint.add_proteins(
                        network,
                        interaction_detection_methods=configuration[
                            "protein-protein interactions"]["MINT"].get(
                                "interaction detection methods", []),
                        interaction_types=configuration[
                            "protein-protein interactions"]["MINT"].get(
                                "interaction types", []),
                        mi_score=configuration["protein-protein interactions"]
                        ["MINT"].get("MI score", 0.0),
                        taxon_identifier=configuration[
                            "protein-protein interactions"]["MINT"].get(
                                "taxon identifier", 9606),
                    )

                if "STRING" in configuration[
                        "protein-protein interactions"] and configuration[
                            "protein-protein interactions"]["STRING"].get(
                                "neighbors", 0) > k:
                    string.add_proteins(
                        network,
                        neighborhood=configuration[
                            "protein-protein interactions"]["STRING"].get(
                                "neighborhood", 0.0),
                        neighborhood_transferred=configuration[
                            "protein-protein interactions"]["STRING"].get(
                                "neighborhood transferred", 0.0),
                        fusion=configuration["protein-protein interactions"]
                        ["STRING"].get("fusion", 0.0),
                        cooccurence=configuration[
                            "protein-protein interactions"]["STRING"].get(
                                "cooccurence", 0.0),
                        homology=configuration["protein-protein interactions"]
                        ["STRING"].get("homology", 0.0),
                        coexpression=configuration[
                            "protein-protein interactions"]["STRING"].get(
                                "coexpression", 0.0),
                        coexpression_transferred=configuration[
                            "protein-protein interactions"]["STRING"].get(
                                "coexpression transferred", 0.0),
                        experiments=configuration[
                            "protein-protein interactions"]["STRING"].get(
                                "experiments", 0.0),
                        experiments_transferred=configuration[
                            "protein-protein interactions"]["STRING"].get(
                                "experiments transferred", 0.0),
                        database=configuration["protein-protein interactions"]
                        ["STRING"].get("database", 0.0),
                        database_transferred=configuration[
                            "protein-protein interactions"]["STRING"].get(
                                "database transferred", 0.0),
                        textmining=configuration[
                            "protein-protein interactions"]["STRING"].get(
                                "textmining", 0.0),
                        textmining_transferred=configuration[
                            "protein-protein interactions"]["STRING"].get(
                                "textmining transferred", 0.0),
                        combined_score=configuration[
                            "protein-protein interactions"]["STRING"].get(
                                "combined score", 0.0),
                        physical=configuration["protein-protein interactions"]
                        ["STRING"].get("physical", False),
                        taxon_identifier=configuration[
                            "protein-protein interactions"]["STRING"].get(
                                "taxon identifier", 9606),
                    )

                k += 1

            if k:
                protein_interaction_network.annotate_proteins(network)
                protein_interaction_network.remove_unannotated_proteins(
                    network)

            if "BioGRID" in configuration["protein-protein interactions"]:
                biogrid.add_interactions(
                    network,
                    interaction_throughput=configuration[
                        "protein-protein interactions"]["BioGRID"].get(
                            "interaction throughput", []),
                    experimental_system=configuration[
                        "protein-protein interactions"]["BioGRID"].get(
                            "experimental system", []),
                    experimental_system_type=configuration[
                        "protein-protein interactions"]["BioGRID"].get(
                            "experimental system type", []),
                    multi_validated_physical=configuration[
                        "protein-protein interactions"]["BioGRID"].get(
                            "multi-validated physical", False),
                    taxon_identifier=configuration[
                        "protein-protein interactions"]["BioGRID"].get(
                            "taxon identifier", 9606),
                )

            if "CORUM" in configuration["protein-protein interactions"]:
                corum.add_interactions(
                    network,
                    protein_complex_purification_method=configuration[
                        "protein-protein interactions"]["CORUM"].get(
                            "protein complex purification method",
                            [],
                        ),
                    taxon_identifier=configuration[
                        "protein-protein interactions"]["CORUM"].get(
                            "taxon identifier", 9606),
                )

            if "IntAct" in configuration["protein-protein interactions"]:
                intact.add_interactions(
                    network,
                    interaction_detection_methods=configuration[
                        "protein-protein interactions"]["IntAct"].get(
                            "interaction detection methods", []),
                    interaction_types=configuration[
                        "protein-protein interactions"]["IntAct"].get(
                            "interaction types", []),
                    mi_score=configuration["protein-protein interactions"]
                    ["IntAct"].get("MI score", 0.0),
                    taxon_identifier=configuration[
                        "protein-protein interactions"]["IntAct"].get(
                            "taxon identifier", 9606),
                )

            if "MINT" in configuration["protein-protein interactions"]:
                mint.add_interactions(
                    network,
                    interaction_detection_methods=configuration[
                        "protein-protein interactions"]["MINT"].get(
                            "interaction detection methods", []),
                    interaction_types=configuration[
                        "protein-protein interactions"]["MINT"].get(
                            "interaction types", []),
                    mi_score=configuration["protein-protein interactions"]
                    ["MINT"].get("MI score", 0.0),
                    taxon_identifier=configuration[
                        "protein-protein interactions"]["MINT"].get(
                            "taxon identifier", 9606),
                )

            if "STRING" in configuration["protein-protein interactions"]:
                string.add_interactions(
                    network,
                    neighborhood=configuration["protein-protein interactions"]
                    ["STRING"].get("neighborhood", 0.0),
                    neighborhood_transferred=configuration[
                        "protein-protein interactions"]["STRING"].get(
                            "neighborhood transferred", 0.0),
                    fusion=configuration["protein-protein interactions"]
                    ["STRING"].get("fusion", 0.0),
                    cooccurence=configuration["protein-protein interactions"]
                    ["STRING"].get("cooccurence", 0.0),
                    homology=configuration["protein-protein interactions"]
                    ["STRING"].get("homology", 0.0),
                    coexpression=configuration["protein-protein interactions"]
                    ["STRING"].get("coexpression", 0.0),
                    coexpression_transferred=configuration[
                        "protein-protein interactions"]["STRING"].get(
                            "coexpression transferred", 0.0),
                    experiments=configuration["protein-protein interactions"]
                    ["STRING"].get("experiments", 0.0),
                    experiments_transferred=configuration[
                        "protein-protein interactions"]["STRING"].get(
                            "experiments transferred", 0.0),
                    database=configuration["protein-protein interactions"]
                    ["STRING"].get("database", 0.0),
                    database_transferred=configuration[
                        "protein-protein interactions"]["STRING"].get(
                            "database transferred", 0.0),
                    textmining=configuration["protein-protein interactions"]
                    ["STRING"].get("textmining", 0.0),
                    textmining_transferred=configuration[
                        "protein-protein interactions"]["STRING"].get(
                            "textmining transferred", 0.0),
                    combined_score=configuration[
                        "protein-protein interactions"]["STRING"].get(
                            "combined score", 0.0),
                    physical=configuration["protein-protein interactions"]
                    ["STRING"].get("physical", False),
                    taxon_identifier=configuration[
                        "protein-protein interactions"]["STRING"].get(
                            "taxon identifier", 9606),
                )

        if "Cytoscape" in configuration:
            if (configuration["Cytoscape"].get("bar chart",
                                               {}).get("type") == "z-score"):
                cytoscape_styles = styles.get_styles(
                    network,
                    bar_chart_range=configuration["Cytoscape"]
                    ["bar chart"].get("range", (-2.0, 2.0)),
                    get_bar_chart_range=protein_interaction_network.
                    get_z_score_range,
                    site_combination=combination.SITE_COMBINATION.get(
                        configuration["Cytoscape"]["bar chart"].get(
                            "combine sites"),
                        combination.SITE_COMBINATION["absmax"],
                    ),
                )

            elif (configuration["Cytoscape"].get(
                    "bar chart", {}).get("type") == "proportion"):
                cytoscape_styles = styles.get_styles(
                    network,
                    bar_chart_range=configuration["Cytoscape"]
                    ["bar chart"].get("range", (0.025, 0.975)),
                    get_bar_chart_range=protein_interaction_network.
                    get_propotion_range,
                    site_combination=combination.SITE_COMBINATION.get(
                        configuration["Cytoscape"]["bar chart"].get(
                            "combine sites"),
                        combination.SITE_COMBINATION["absmax"],
                    ),
                )

            else:
                cytoscape_styles = styles.get_styles(
                    network,
                    bar_chart_range=configuration["Cytoscape"]
                    ["bar chart"].get("range", (-1.0, 1.0)),
                    site_combination=combination.SITE_COMBINATION.get(
                        configuration["Cytoscape"]["bar chart"].get(
                            "combine sites"),
                        combination.SITE_COMBINATION["absmax"],
                    ),
                )

            styles.export(
                cytoscape_styles,
                os.path.splitext(os.path.basename(configuration_file))[0],
                ".{}".format(i) if len(configurations) > 1 else "",
            )

            protein_interaction_network.set_post_translational_modification(
                network)

            if (configuration["Cytoscape"].get("node color",
                                               {}).get("type") == "z-score"):
                protein_interaction_network.set_changes(
                    network,
                    site_combination=combination.SITE_COMBINATION[
                        configuration["Cytoscape"]["node color"].get(
                            "combine sites", "absmax")],
                    changes=configuration["Cytoscape"]["node color"].get(
                        "range", (-2.0, 2.0)),
                    get_range=protein_interaction_network.get_z_score_range,
                )

            elif (configuration["Cytoscape"].get(
                    "node color", {}).get("type") == "proportion"):
                protein_interaction_network.set_changes(
                    network,
                    site_combination=combination.SITE_COMBINATION[
                        configuration["Cytoscape"]["node color"].get(
                            "combine sites", "absmax")],
                    changes=configuration["Cytoscape"]["node color"].get(
                        "range", (0.025, 0.975)),
                    get_range=protein_interaction_network.get_proportion_range,
                )

            else:
                protein_interaction_network.set_changes(
                    network,
                    site_combination=combination.SITE_COMBINATION[
                        configuration["Cytoscape"]["node color"].get(
                            "combine sites", "absmax")],
                    changes=configuration["Cytoscape"]["node color"].get(
                        "range", (-1.0, 1.0)),
                )

        protein_interaction_network.export(
            network,
            os.path.splitext(os.path.basename(configuration_file))[0],
            ".{}".format(i) if len(configurations) > 1 else "",
        )

        if "post-processing" in configuration:
            if "neighborhood extraction" in configuration:
                for protein in configuration["neighborhood"].get(
                        "proteins", []):
                    protein_interaction_network.export(
                        protein_interaction_network.get_neighborhood(
                            network,
                            protein,
                            configuration["neighborhood"].get("distance", 1),
                            configuration["neighborhood"].get(
                                "isoforms", True),
                        ),
                        os.path.splitext(
                            os.path.basename(configuration_file))[0],
                        ".{}.{}.{}".format(i, protein,
                                           network.nodes[protein]["gene name"])
                        if len(configurations) > 1 else ".{}.{}".format(
                            protein, network.nodes[protein]["gene name"]),
                    )

            if "module detection" in configuration["post-processing"]:
                protein_interaction_network.set_edge_weights(
                    network,
                    weight=combination.CONFIDENCE_SCORE_COMBINATION[
                        configuration["post-processing"]
                        ["module detection"].get("edge weight", "number")])

                for j, module in enumerate(
                        protein_interaction_network.get_modules(
                            network,
                            module_size=configuration["post-processing"]
                            ["module detection"].get(
                                "module size", network.number_of_nodes()),
                            module_size_combination=combination.
                            MODULE_SIZE_COMBINATION.get(
                                configuration["post-processing"]
                                ["module detection"].get("combine sizes"),
                                combination.MODULE_SIZE_COMBINATION["mean"],
                            ),
                            algorithm=algorithm.ALGORITHM.get(
                                configuration["post-processing"]
                                ["module detection"].get("algorithm"),
                                algorithm.ALGORITHM["Louvain"],
                            ),
                        ),
                        start=1,
                ):
                    protein_interaction_network.export(
                        network.subgraph(module),
                        os.path.splitext(
                            os.path.basename(configuration_file))[0],
                        ".{}.{}".format(i, j)
                        if len(configurations) > 1 else ".{}".format(j),
                    )

            if "enrichment analysis" in configuration["post-processing"]:
                protein_interaction_network.set_edge_weights(
                    network,
                    weight=combination.CONFIDENCE_SCORE_COMBINATION[
                        configuration["post-processing"]
                        ["enrichment analysis"].get("edge weight", "number")])

                if (configuration["post-processing"]
                    ["enrichment analysis"].get("type") == "z-score"):
                    modules, p_values = protein_interaction_network.get_change_enriched_modules(
                        network,
                        module_size=configuration["post-processing"]
                        ["enrichment analysis"].get("module size",
                                                    network.number_of_nodes()),
                        p=configuration["post-processing"]
                        ["enrichment analysis"].get("p", 0.05),
                        changes=configuration["post-processing"]
                        ["enrichment analysis"].get("range", (-2.0, 2.0)),
                        get_range=protein_interaction_network.
                        get_z_score_range,
                        site_combination=combination.SITE_COMBINATION.get(
                            configuration["post-processing"]
                            ["enrichment analysis"].get("combine sites"),
                            combination.SITE_COMBINATION["absmax"],
                        ),
                        module_size_combination=combination.
                        MODULE_SIZE_COMBINATION.get(
                            configuration["post-processing"]
                            ["enrichment analysis"].get("combine sites"),
                            combination.MODULE_SIZE_COMBINATION["mean"],
                        ),
                        test=configuration["post-processing"]
                        ["enrichment analysis"].get("test", "outside"),
                        algorithm=algorithm.ALGORITHM.get(
                            configuration["post-processing"]
                            ["enrichment analysis"].get("algorithm"),
                            algorithm.ALGORITHM["Louvain"],
                        ),
                    )

                elif (configuration["post-processing"]
                      ["enrichment analysis"].get("type") == "proportion"):
                    modules, p_values = protein_interaction_network.get_change_enriched_modules(
                        network,
                        module_size=configuration["post-processing"]
                        ["enrichment analysis"].get("module size",
                                                    network.number_of_nodes()),
                        p=configuration["post-processing"]
                        ["enrichment analysis"].get("p", 0.05),
                        changes=configuration["post-processing"]
                        ["enrichment analysis"].get("range", (0.025, 0.975)),
                        get_range=protein_interaction_network.
                        get_proportion_range,
                        site_combination=combination.SITE_COMBINATION.get(
                            configuration["post-processing"]
                            ["enrichment analysis"].get("combine sites"),
                            combination.SITE_COMBINATION["absmax"],
                        ),
                        module_size_combination=combination.
                        MODULE_SIZE_COMBINATION.get(
                            configuration["post-processing"]
                            ["enrichment analysis"].get("combine sites"),
                            combination.MODULE_SIZE_COMBINATION["mean"],
                        ),
                        test=configuration["post-processing"]
                        ["enrichment analysis"].get("test", "outside"),
                        algorithm=algorithm.ALGORITHM.get(
                            configuration["post-processing"]
                            ["enrichment analysis"].get("algorithm"),
                            algorithm.ALGORITHM["Louvain"],
                        ),
                    )

                else:
                    modules, p_values = protein_interaction_network.get_change_enriched_modules(
                        network,
                        module_size=configuration["post-processing"]
                        ["enrichment analysis"].get("module size",
                                                    network.number_of_nodes()),
                        p=configuration["post-processing"]
                        ["enrichment analysis"].get("p", 0.05),
                        changes=configuration["post-processing"]
                        ["enrichment analysis"].get("range", (-1.0, 1.0)),
                        site_combination=combination.SITE_COMBINATION.get(
                            configuration["post-processing"]
                            ["enrichment analysis"].get("combine sites"),
                            combination.SITE_COMBINATION["absmax"],
                        ),
                        module_size_combination=combination.
                        MODULE_SIZE_COMBINATION.get(
                            configuration["post-processing"]
                            ["enrichment analysis"].get("combine sites"),
                            combination.MODULE_SIZE_COMBINATION["mean"],
                        ),
                        test=configuration["post-processing"]
                        ["enrichment analysis"].get("test", "outside"),
                        algorithm=algorithm.ALGORITHM.get(
                            configuration["post-processing"]
                            ["enrichment analysis"].get("algorithm"),
                            algorithm.ALGORITHM["Louvain"],
                        ),
                    )

                protein_interaction_network.remove_edge_weights(network)

                for time in p_values:
                    for modification in sorted(p_values[time]):
                        for j in sorted(p_values[time][modification]):
                            logger.info("{}\t{}\t{}\t{:.2e}".format(
                                time,
                                modification,
                                j + 1,
                                p_values[time][modification][j],
                            ))
                            protein_interaction_network.export(
                                network.subgraph(modules[j]),
                                os.path.splitext(
                                    os.path.basename(configuration_file))[0],
                                ".{}.{}".format(i, j + 1) if
                                len(configurations) > 1 else ".{}".format(j +
                                                                          1),
                            )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-c",
        "--configurations",
        help="configuration files",
        nargs="+",
        required=True,
    )
    parser.add_argument(
        "-l",
        "--log",
        action="store_true",
        default=False,
        help="store a log file for each configuration file (default: False)")

    parser.add_argument(
        "-p",
        "--processes",
        help="maximum number of processes (default: {})".format(
            os.cpu_count()),
        type=int,
        default=os.cpu_count())
    args = parser.parse_args()

    with concurrent.futures.ProcessPoolExecutor(
            max_workers=args.processes) as executor:
        executor.map(process, args.configurations, itertools.repeat(args.log))


if __name__ == "__main__":
    main()
