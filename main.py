import argparse
import json
import logging
import os
import sys
import xml.etree.ElementTree as ET

import networkx as nx

from pipeline.cytoscape_styles import CytoscapeStyles
from pipeline.interface import convert, merge, extract
from pipeline.protein_protein_interaction_network import (
    ProteinProteinInteractionNetwork,
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
            network = ProteinProteinInteractionNetwork()

            for entry in configuration.get("PTM", {}):
                for (
                    gene_name,
                    protein,
                    protein_name,
                    sites,
                ) in network.add_proteins_from_spreadsheet(
                    entry["file"],
                    entry["label"],
                    entry["time"],
                    protein_accession_column=entry["protein column"],
                    position_col=entry["position column"],
                    replicates=entry["replicate columns"],
                    protein_accession_format=extract.EXTRACT.get(
                        entry["protein format"], lambda x: x
                    ),
                    position_format=extract.EXTRACT.get(
                        entry["position format"], lambda x: x
                    ),
                    sheet_name=entry.get("sheet", 0),
                    header=entry.get("header", 1) - 1,
                    num_sites=entry.get("sites", 5),
                    num_replicates=entry.get("replicates", 2),
                    merge_replicates=merge.MERGE.get(
                        entry.get("merge replicates", "mean"), merge.MERGE["mean"]
                    ),
                    convert_measurement=convert.LOG_BASE[entry.get("log base")],
                ):
                    logger.info(
                        "{}\t{}\t{}\t{}\t{}\t{}".format(
                            gene_name,
                            protein_name,
                            protein if "-" in protein else "{}  ".format(protein),
                            entry["time"],
                            entry["label"],
                            " ".join(["{:.2f}".format(site) for site in sites]),
                        )
                    )

            if configuration.get("PPI"):
                if "BioGRID" in configuration["PPI"]:
                    for (
                        interactor_a,
                        interactor_b,
                        score,
                    ) in network.add_interactions_from_biogrid(
                        experimental_system=configuration["PPI"]["BioGRID"].get(
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
                        multi_validated_physical=configuration["PPI"]["BioGRID"].get(
                            "multi-validated physical", False
                        ),
                    ):
                        logger.info(
                            "{}\t{}\tBioGRID\t{:.3f}".format(
                                *sorted(
                                    [
                                        interactor_a
                                        if "-" in interactor_a
                                        else "{}  ".format(interactor_a),
                                        interactor_b
                                        if "-" in interactor_b
                                        else "{}  ".format(interactor_b),
                                    ]
                                ),
                                score
                            )
                        )

                if "CORUM" in configuration["PPI"]:
                    for (
                        interactor_a,
                        interactor_b,
                        score,
                    ) in network.add_interactions_from_corum(
                        protein_complex_purification_method=configuration["PPI"][
                            "CORUM"
                        ].get(
                            "protein complex purification method",
                            [],
                        )
                    ):
                        logger.info(
                            "{}\t{}\tCORUM\t{:.3f}".format(
                                *sorted(
                                    [
                                        interactor_a
                                        if "-" in interactor_a
                                        else "{}  ".format(interactor_a),
                                        interactor_b
                                        if "-" in interactor_b
                                        else "{}  ".format(interactor_b),
                                    ]
                                ),
                                score
                            )
                        )

                if "IntAct" in configuration["PPI"]:
                    for (
                        interactor_a,
                        interactor_b,
                        score,
                    ) in network.add_interactions_from_intact(
                        interaction_detection_methods=configuration["PPI"][
                            "IntAct"
                        ].get("interaction detection methods", []),
                        interaction_types=configuration["PPI"]["IntAct"].get(
                            "interaction_types", []
                        ),
                        mi_score=configuration["PPI"]["IntAct"].get("MI score", 0.27),
                    ):
                        logger.info(
                            "{}\t{}\tIntAct\t{:.3f}".format(
                                *sorted(
                                    [
                                        interactor_a
                                        if "-" in interactor_a
                                        else "{}  ".format(interactor_a),
                                        interactor_b
                                        if "-" in interactor_b
                                        else "{}  ".format(interactor_b),
                                    ]
                                ),
                                score
                            )
                        )

                if "Reactome" in configuration["PPI"]:
                    for (
                        interactor_a,
                        interactor_b,
                        score,
                    ) in network.add_interactions_from_reactome(
                        interaction_detection_methods=configuration["PPI"][
                            "Reactome"
                        ].get("interaction detection methods", []),
                        interaction_types=configuration["PPI"]["Reactome"].get(
                            "interaction_types", []
                        ),
                        reactome_score=configuration["PPI"]["Reactome"].get(
                            "Reactome score", 0.0
                        ),
                    ):
                        logger.info(
                            "{}\t{}\tReactome\t{:.3f}".format(
                                *sorted(
                                    [
                                        interactor_a
                                        if "-" in interactor_a
                                        else "{}  ".format(interactor_a),
                                        interactor_b
                                        if "-" in interactor_b
                                        else "{}  ".format(interactor_b),
                                    ]
                                ),
                                score
                            )
                        )

                if "STRING" in configuration["PPI"]:
                    for (
                        interactor_a,
                        interactor_b,
                        score,
                    ) in network.add_interactions_from_string(
                        neighborhood=configuration["PPI"]["STRING"].get(
                            "neighborhood", 0.0
                        ),
                        neighborhood_transferred=configuration["PPI"]["STRING"].get(
                            "neighborhood transferred", 0.0
                        ),
                        fusion=configuration["PPI"]["STRING"].get("fusion", 0.0),
                        cooccurence=configuration["PPI"]["STRING"].get(
                            "cooccurence", 0.0
                        ),
                        homology=configuration["PPI"]["STRING"].get("homology", 0.0),
                        coexpression=configuration["PPI"]["STRING"].get(
                            "coexpression", 0.0
                        ),
                        coexpression_transferred=configuration["PPI"]["STRING"].get(
                            "coexpression transferred", 0.0
                        ),
                        experiments=configuration["PPI"]["STRING"].get(
                            "experiments", 0.7
                        ),
                        experiments_transferred=configuration["PPI"]["STRING"].get(
                            "experiments transferred", 0.0
                        ),
                        database=configuration["PPI"]["STRING"].get("database", 0.0),
                        database_transferred=configuration["PPI"]["STRING"].get(
                            "database transferred", 0.0
                        ),
                        textmining=configuration["PPI"]["STRING"].get(
                            "textmining", 0.0
                        ),
                        textmining_transferred=configuration["PPI"]["STRING"].get(
                            "textmining transferred", 0.0
                        ),
                        combined_score=configuration["PPI"]["STRING"].get(
                            "combined score", 0.7
                        ),
                        physical=configuration["PPI"]["STRING"].get("physical", False),
                    ):
                        logger.info(
                            "{}\t{}\tSTRING\t{:.3f}".format(
                                *sorted(
                                    [
                                        interactor_a
                                        if "-" in interactor_a
                                        else "{}  ".format(interactor_a),
                                        interactor_b
                                        if "-" in interactor_b
                                        else "{}  ".format(interactor_b),
                                    ]
                                ),
                                score
                            )
                        )

            network.set_post_translational_modification()

            if configuration.get("Cytoscape", {}).get("type") == "z-score":
                network.set_change_data(
                    merge_sites=merge.MERGE[
                        configuration.get("Cytoscape", {}).get("merge sites", "mean")
                    ],
                    change=configuration["Cytoscape"].get("threshold", 2.0),
                    get_range=network.get_z_score_range,
                )

            elif configuration.get("Cytoscape", {}).get("type") == "p-value":
                network.set_change_data(
                    merge_sites=merge.MERGE[
                        configuration.get("Cytoscape", {}).get("merge sites", "mean")
                    ],
                    change=configuration["Cytoscape"].get("threshold", 0.05),
                    get_range=network.get_p_value_range,
                )

            else:
                network.set_change_data(
                    merge_sites=merge.MERGE[
                        configuration.get("Cytoscape", {}).get("merge sites", "mean")
                    ],
                    change=configuration["Cytoscape"].get("threshold", 1.0),
                    get_range=lambda time, ptm, change, merge_sites: (-change, change),
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
                    bar_chart_range=configuration.get("Cytoscape", {})
                    .get("bar chart", {})
                    .get("range", [-3.0, 3.0]),
                    get_bar_chart_range=network.get_z_score_range,
                    merge_sites=merge.MERGE.get(
                        configuration.get("Cytoscape", {}).get("merge sites"),
                        merge.MERGE["mean"],
                    ),
                )
            elif (
                configuration.get("Cytoscape", {}).get("bar chart", {}).get("type")
                == "p-value"
            ):
                styles = CytoscapeStyles(
                    network,
                    bar_chart_range=configuration.get("Cytoscape", {})
                    .get("bar chart", {})
                    .get("range", [5, 95]),
                    get_bar_chart_range=network.get_p_value_range,
                    merge_sites=merge.MERGE.get(
                        configuration.get("Cytoscape", {}).get("merge sites"),
                        merge.MERGE["mean"],
                    ),
                )
            else:
                styles = CytoscapeStyles(
                    network,
                    bar_chart_range=configuration.get("Cytoscape", {})
                    .get("bar chart", {})
                    .get("range", [-3.0, 3.0]),
                    get_bar_chart_range=lambda time, ptm, change, merge_sites: (
                        -change,
                        change,
                    ),
                    merge_sites=merge.MERGE.get(
                        configuration.get("Cytoscape", {}).get("merge sites"),
                        merge.MERGE["mean"],
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

            if configuration.get("module change enrichment"):
                network.set_edge_weights(
                    weight=lambda confidence_scores: len(confidence_scores)
                )
                if configuration["module change enrichment"].get("type") == "z-score":
                    modules, p_values = network.get_module_change_enrichment(
                        p=configuration["module change enrichment"].get("p", 0.05),
                        change=configuration["module change enrichment"].get(
                            "threshold", 2.0
                        ),
                        get_range=network.get_z_score_range,
                        merge_sites=merge.MERGE.get(
                            configuration["module change enrichment"].get(
                                "merge sites"
                            ),
                            merge.MERGE["mean"],
                        ),
                        module_size=configuration["module change enrichment"].get(
                            "module size", (3, 100)
                        ),
                    )

                elif configuration["module change enrichment"].get("type") == "p-value":
                    modules, p_values = network.get_module_change_enrichment(
                        p=configuration["module change enrichment"].get("p", 0.05),
                        change=configuration["module change enrichment"].get(
                            "threshold", 0.05
                        ),
                        get_range=network.get_p_value_range,
                        merge_sites=merge.MERGE.get(
                            configuration["module change enrichment"].get(
                                "merge sites"
                            ),
                            merge.MERGE["mean"],
                        ),
                        module_size=configuration["module change enrichment"].get(
                            "module size", (3, 100)
                        ),
                    )

                else:
                    modules, p_values = network.get_module_change_enrichment(
                        p=configuration["module change enrichment"].get("p", 0.05),
                        change=configuration["module change enrichment"].get(
                            "threshold", 1.0
                        ),
                        get_range=lambda time, ptm, change, merge_sites: (
                            -change,
                            change,
                        ),
                        merge_sites=merge.MERGE.get(
                            configuration["module change enrichment"].get(
                                "merge sites"
                            ),
                            merge.MERGE["mean"],
                        ),
                        module_size=configuration["module change enrichment"].get(
                            "module size", (3, 100)
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
                                os.path.splitext(os.path.basename(configuration_file))[
                                    0
                                ],
                                ".{}.{}".format(i, module + 1)
                                if len(configurations) > 1
                                else ".{}".format(module + 1),
                            )


if __name__ == "__main__":
    main()
