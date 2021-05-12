import argparse
import json
import logging
import os
import sys
import xml.etree.ElementTree as ET

import networkx as nx

from pipeline.cytoscape_styles import CytoscapeStyles
from pipeline.interface import algorithm, convert, combine, extract
from pipeline.protein_protein_interaction_network import (
    ProteinInteractionNetwork,
)


def export_network(network, basename, suffix=""):
    nx.write_graphml_xml(network, "{0}{1}.graphml".format(basename, suffix))


def export_styles(styles, basename, suffix=""):
    styles.write(
        "{0}{1}.xml".format(basename, suffix),
        encoding="UTF-8",
        xml_declaration=True,
    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-c",
        "--configurations",
        help="JSON configuration files specifying protein-protein interaction network assembly and analysis",
        nargs="+",
        required=True,
    )
    parser.add_argument("-l", "--log", help="file name for log")
    args = parser.parse_args()

    logging.basicConfig(stream=sys.stdout, level=logging.DEBUG, format="%(message)s")

    logger = logging.getLogger("root")

    if args.log:
        log_file = logging.FileHandler(args.log, mode="w")
        log_file.setFormatter(logging.Formatter(fmt="%(message)s"))
        logger.addHandler(log_file)

    for configuration_file in args.configurations:
        with open(configuration_file) as configuration:
            configurations = json.load(configuration)

        for i, configuration in enumerate(configurations, start=1):
            network = ProteinInteractionNetwork()

            for entry in configuration.get("genes", {}):
                for (
                    gene_name,
                    protein,
                    protein_name,
                ) in network.add_genes_from_spreadsheet(
                    entry["file"],
                    gene_accession_column=entry["accession column"],
                    gene_accession_format=extract.EXTRACT.get(
                        entry.get("accession format"), lambda x: [x]
                    ),
                    sheet_name=entry.get("sheet", 0),
                    header=entry.get("header", 1) - 1,
                    organism=entry.get("organism", 9606),
                ):
                    logger.info("{}\t{}\t{}".format(protein, gene_name, protein_name))

            for entry in configuration.get("proteins", {}):
                for (
                    gene_name,
                    protein,
                    protein_name,
                ) in network.add_proteins_from_spreadsheet(
                    entry["file"],
                    protein_accession_column=entry["accession column"],
                    protein_accession_format=extract.EXTRACT.get(
                        entry.get("accession format"), lambda x: [x]
                    ),
                    time=entry.get("time", 0),
                    ptm=entry.get("post-translational modification", ""),
                    position_column=entry.get("position column", ""),
                    position_format=extract.EXTRACT.get(
                        entry.get("position format"), lambda x: [x]
                    ),
                    replicates=entry.get("replicate columns", []),
                    sheet_name=entry.get("sheet", 0),
                    header=entry.get("header", 1) - 1,
                    num_sites=entry.get("sites", 1000),
                    num_replicates=entry.get("replicates", 1),
                    combine_replicates=combine.COMBINE.get(
                        entry.get("combine replicates", "mean"), combine.COMBINE["mean"]
                    ),
                    convert_measurement=convert.LOG_BASE[entry.get("log base")],
                    organism=entry.get("organism", 9606),
                ):
                    logger.info("{}\t{}\t{}".format(protein, gene_name, protein_name))

            if "protein-protein interactions" in configuration:
                if "BioGRID" in configuration["protein-protein interactions"]:
                    for (
                        interactor_a,
                        interactor_b,
                        score,
                    ) in network.add_interactions_from_biogrid(
                        experimental_system=configuration[
                            "protein-protein interactions"
                        ]["BioGRID"].get(
                            "experimental system",
                            [
                                "Affinity Capture-Luminescence",
                                "Affinity Capture-MS",
                                "Affinity Capture-RNA",
                                "Affinity Capture-Western",
                                "Biochemical Activity",
                                "Co-crystal Structure",
                                "Co-purification",
                                "FRET",
                                "PCA",
                                "Two-hybrid",
                            ],
                        ),
                        organism=configuration["protein-protein interactions"][
                            "BioGRID"
                        ].get("organism", 9606),
                        multi_validated_physical=configuration[
                            "protein-protein interactions"
                        ]["BioGRID"].get("multi-validated physical", False),
                    ):
                        logger.info(
                            "{}\t{}\tBioGRID\t{:.3f}".format(
                                *sorted(
                                    [
                                        interactor_a,
                                        interactor_b,
                                    ]
                                ),
                                score
                            )
                        )

                if "CORUM" in configuration["protein-protein interactions"]:
                    for (
                        interactor_a,
                        interactor_b,
                        score,
                    ) in network.add_interactions_from_corum(
                        protein_complex_purification_method=configuration[
                            "protein-protein interactions"
                        ]["CORUM"].get(
                            "protein complex purification method",
                            [],
                        ),
                        organism=configuration["protein-protein interactions"][
                            "CORUM"
                        ].get("organism", "Human"),
                    ):
                        logger.info(
                            "{}\t{}\tCORUM\t{:.3f}".format(
                                *sorted([interactor_a, interactor_b]), score
                            )
                        )

                if "IntAct" in configuration["protein-protein interactions"]:
                    for (
                        interactor_a,
                        interactor_b,
                        score,
                    ) in network.add_interactions_from_intact(
                        interaction_detection_methods=configuration[
                            "protein-protein interactions"
                        ]["IntAct"].get("interaction detection methods", []),
                        interaction_types=configuration["protein-protein interactions"][
                            "IntAct"
                        ].get("interaction_types", []),
                        mi_score=configuration["protein-protein interactions"][
                            "IntAct"
                        ].get("MI score", 0.27),
                        organism=configuration["protein-protein interactions"][
                            "IntAct"
                        ].get("organism", 9606),
                    ):
                        logger.info(
                            "{}\t{}\tIntAct\t{:.3f}".format(
                                *sorted([interactor_a, interactor_b]), score
                            )
                        )

                if "Reactome" in configuration["protein-protein interactions"]:
                    for (
                        interactor_a,
                        interactor_b,
                        score,
                    ) in network.add_interactions_from_reactome(
                        interaction_detection_methods=configuration[
                            "protein-protein interactions"
                        ]["Reactome"].get("interaction detection methods", []),
                        interaction_types=configuration["protein-protein interactions"][
                            "Reactome"
                        ].get("interaction_types", []),
                        reactome_score=configuration["protein-protein interactions"][
                            "Reactome"
                        ].get("Reactome score", 0.0),
                        organism=configuration["protein-protein interactions"][
                            "Reactome"
                        ].get("organism", 9606),
                    ):
                        logger.info(
                            "{}\t{}\tReactome\t{:.3f}".format(
                                *sorted([interactor_a, interactor_b]), score
                            )
                        )

                if "STRING" in configuration["protein-protein interactions"]:
                    for (
                        interactor_a,
                        interactor_b,
                        score,
                    ) in network.add_interactions_from_string(
                        neighborhood=configuration["protein-protein interactions"][
                            "STRING"
                        ].get("neighborhood", 0.0),
                        neighborhood_transferred=configuration[
                            "protein-protein interactions"
                        ]["STRING"].get("neighborhood transferred", 0.0),
                        fusion=configuration["protein-protein interactions"][
                            "STRING"
                        ].get("fusion", 0.0),
                        cooccurence=configuration["protein-protein interactions"][
                            "STRING"
                        ].get("cooccurence", 0.0),
                        homology=configuration["protein-protein interactions"][
                            "STRING"
                        ].get("homology", 0.0),
                        coexpression=configuration["protein-protein interactions"][
                            "STRING"
                        ].get("coexpression", 0.0),
                        coexpression_transferred=configuration[
                            "protein-protein interactions"
                        ]["STRING"].get("coexpression transferred", 0.0),
                        experiments=configuration["protein-protein interactions"][
                            "STRING"
                        ].get("experiments", 0.7),
                        experiments_transferred=configuration[
                            "protein-protein interactions"
                        ]["STRING"].get("experiments transferred", 0.0),
                        database=configuration["protein-protein interactions"][
                            "STRING"
                        ].get("database", 0.0),
                        database_transferred=configuration[
                            "protein-protein interactions"
                        ]["STRING"].get("database transferred", 0.0),
                        textmining=configuration["protein-protein interactions"][
                            "STRING"
                        ].get("textmining", 0.0),
                        textmining_transferred=configuration[
                            "protein-protein interactions"
                        ]["STRING"].get("textmining transferred", 0.0),
                        combined_score=configuration["protein-protein interactions"][
                            "STRING"
                        ].get("combined score", 0.7),
                        organism=configuration["protein-protein interactions"][
                            "STRING"
                        ].get("organism", 9606),
                        physical=configuration["protein-protein interactions"][
                            "STRING"
                        ].get("physical", False),
                    ):
                        logger.info(
                            "{}\t{}\tSTRING\t{:.3f}".format(
                                *sorted([interactor_a, interactor_b]), score
                            )
                        )

            if "Cytoscape" in configuration:
                network.set_post_translational_modification()

                if (
                    configuration["Cytoscape"].get("node color", {}).get("type")
                    == "z-score"
                ):
                    network.set_changes(
                        combine_sites=combine.COMBINE[
                            configuration["Cytoscape"]["node color"].get(
                                "combine sites", "mean"
                            )
                        ],
                        changes=configuration["Cytoscape"]["node color"].get(
                            "range", (-2.0, 2.0)
                        ),
                        get_range=network.get_z_score_range,
                    )

                elif (
                    configuration["Cytoscape"].get("node color", {}).get("type")
                    == "proportion"
                ):
                    network.set_changes(
                        combine_sites=combine.COMBINE[
                            configuration["Cytoscape"]["node color"].get(
                                "combine sites", "mean"
                            )
                        ],
                        changes=configuration["Cytoscape"]["node color"].get(
                            "range", (0.025, 0.975)
                        ),
                        get_range=network.get_proportion_range,
                    )

                else:
                    network.set_changes(
                        combine_sites=combine.COMBINE[
                            configuration["Cytoscape"]["node color"].get(
                                "combine sites", "mean"
                            )
                        ],
                        changes=configuration["Cytoscape"]["node color"].get(
                            "range", (-1.0, 1.0)
                        ),
                    )

            export_network(
                network,
                os.path.splitext(os.path.basename(configuration_file))[0],
                ".{}".format(i) if len(configurations) > 1 else "",
            )

            if (
                configuration.get("Cytoscape", {}).get("bar chart", {}).get("type")
                == "z-score"
            ):
                styles = CytoscapeStyles(
                    network,
                    bar_chart_range=configuration["Cytoscape"]["bar chart"].get(
                        "range", (-2.0, 2.0)
                    ),
                    get_bar_chart_range=network.get_z_score_range,
                    combine_sites=combine.COMBINE.get(
                        configuration["Cytoscape"]["bar chart"].get("combine sites"),
                        combine.COMBINE["mean"],
                    ),
                )

            elif (
                configuration.get("Cytoscape", {}).get("bar chart", {}).get("type")
                == "proportion"
            ):
                styles = CytoscapeStyles(
                    network,
                    bar_chart_range=configuration["Cytoscape"]["bar chart"].get(
                        "range", (0.025, 0.975)
                    ),
                    get_bar_chart_range=network.get_propotion_range,
                    combine_sites=combine.COMBINE.get(
                        configuration["Cytoscape"]["bar chart"].get("combine sites"),
                        combine.COMBINE["mean"],
                    ),
                )

            else:
                styles = CytoscapeStyles(
                    network,
                    bar_chart_range=configuration["Cytoscape"]["bar chart"].get(
                        "range", (-1.0, 1.0)
                    ),
                    combine_sites=combine.COMBINE.get(
                        configuration["Cytoscape"]["bar chart"].get("combine sites"),
                        combine.COMBINE["mean"],
                    ),
                )

            export_styles(
                styles,
                os.path.splitext(os.path.basename(configuration_file))[0],
                ".{}".format(i) if len(configurations) > 1 else "",
            )

            if configuration.get("neighborhood"):
                for protein in configuration["neighborhood"].get("proteins", []):
                    export_network(
                        network.get_neighborhood(
                            protein,
                            configuration["neighborhood"].get("distance", 1),
                            configuration["neighborhood"].get("isoforms", True),
                        ),
                        os.path.splitext(os.path.basename(configuration_file))[0],
                        ".{}.{}.{}".format(
                            i, protein, network.nodes[protein]["gene name"]
                        )
                        if len(configurations) > 1
                        else ".{}.{}".format(
                            protein, network.nodes[protein]["gene name"]
                        ),
                    )

            if "post-processing" in configuration:
                if "enrichment analysis" in configuration["post-processing"]:
                    network.set_edge_weights(
                        weight=lambda confidence_scores: len(confidence_scores)
                    )
                    if (
                        configuration["post-processing"]["enrichment analysis"].get(
                            "type"
                        )
                        == "z-score"
                    ):
                        modules, p_values = network.get_module_change_enrichment(
                            p=configuration["post-processing"][
                                "enrichment analysis"
                            ].get("p", 0.05),
                            changes=configuration["post-processing"][
                                "enrichment analysis"
                            ].get("range", (-2.0, 2.0)),
                            get_range=network.get_z_score_range,
                            combine_sites=combine.COMBINE.get(
                                configuration["post-processing"][
                                    "enrichment analysis"
                                ].get("combine sites"),
                                combine.COMBINE["mean"],
                            ),
                            module_size=configuration["post-processing"][
                                "enrichment analysis"
                            ].get("module size", 35),
                            combine_sizes=combine.COMBINE.get(
                                configuration["post-processing"][
                                    "enrichment analysis"
                                ].get("combine sites"),
                                combine.COMBINE["mean"],
                            ),
                            test=configuration["post-processing"][
                                "enrichment analysis"
                            ].get("test", "two-sided"),
                            algorithm=algorithm.ALGORITHM.get(
                                configuration["post-processing"][
                                    "enrichment analysis"
                                ].get("algorithm"),
                                algorithm.ALGORITHM["Louvain"],
                            ),
                        )

                    elif (
                        configuration["post-processing"]["enrichment analysis"].get(
                            "type"
                        )
                        == "proportion"
                    ):
                        modules, p_values = network.get_module_change_enrichment(
                            p=configuration["post-processing"][
                                "enrichment analysis"
                            ].get("p", 0.05),
                            changes=configuration["post-processing"][
                                "enrichment analysis"
                            ].get("range", (0.025, 0.975)),
                            get_range=network.get_proportion_range,
                            combine_sites=combine.COMBINE.get(
                                configuration["post-processing"][
                                    "enrichment analysis"
                                ].get("combine sites"),
                                combine.COMBINE["mean"],
                            ),
                            module_size=configuration["post-processing"][
                                "enrichment analysis"
                            ].get("module size", 35),
                            combine_sizes=combine.COMBINE.get(
                                configuration["post-processing"][
                                    "enrichment analysis"
                                ].get("combine sites"),
                                combine.COMBINE["mean"],
                            ),
                            test=configuration["post-processing"][
                                "enrichment analysis"
                            ].get("test", "two-sided"),
                            algorithm=algorithm.ALGORITHM.get(
                                configuration["post-processing"][
                                    "enrichment analysis"
                                ].get("algorithm"),
                                algorithm.ALGORITHM["Louvain"],
                            ),
                        )

                    else:
                        modules, p_values = network.get_module_change_enrichment(
                            p=configuration["post-processing"][
                                "enrichment analysis"
                            ].get("p", 0.05),
                            changes=configuration["post-processing"][
                                "enrichment analysis"
                            ].get("range", (-1.0, 1.0)),
                            combine_sites=combine.COMBINE.get(
                                configuration["post-processing"][
                                    "enrichment analysis"
                                ].get("combine sites"),
                                combine.COMBINE["mean"],
                            ),
                            module_size=configuration["post-processing"][
                                "enrichment analysis"
                            ].get("module size", 35),
                            combine_sizes=combine.COMBINE.get(
                                configuration["post-processing"][
                                    "enrichment analysis"
                                ].get("combine sites"),
                                combine.COMBINE["mean"],
                            ),
                            test=configuration["post-processing"][
                                "enrichment analysis"
                            ].get("test", "two-sided"),
                            algorithm=algorithm.ALGORITHM.get(
                                configuration["post-processing"][
                                    "enrichment analysis"
                                ].get("algorithm"),
                                algorithm.ALGORITHM["Louvain"],
                            ),
                        )

                    for time in p_values:
                        for ptm in sorted(p_values[time]):
                            for module in sorted(p_values[time][ptm]):
                                logger.info(
                                    "{}\t{}\t{}\t{:.2E}".format(
                                        time,
                                        ptm,
                                        module + 1,
                                        p_values[time][ptm][module],
                                    )
                                )

                                export_network(
                                    network.subgraph(modules[module]),
                                    os.path.splitext(
                                        os.path.basename(configuration_file)
                                    )[0],
                                    ".{}.{}".format(i, module + 1)
                                    if len(configurations) > 1
                                    else ".{}".format(module + 1),
                                )


if __name__ == "__main__":
    main()
