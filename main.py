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


def export_network(path, network, suffix=""):
    if path.endswith(".graphml") or path.endswith(".xml"):
        nx.write_graphml_xml(
            network, "{0}{2}{1}".format(*os.path.splitext(path), suffix)
        )
    elif path.endswith(".cyjs") or path.endswith(".json"):
        with open("{0}{2}{1}".format(*os.path.splitext(path), suffix), "w") as file:
            json.dump(nx.readwrite.json_graph.cytoscape_data(network), file, indent=2)


def export_styles(path, styles):
    if path.endswith(".xml"):
        with open(path, "wb") as file:
            styles.write(file, encoding="UTF-8", xml_declaration=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-c",
        "--configurations",
        help="YAML configuration files",
        nargs="*",
        required=True,
    )
    parser.add_argument("-l", "--log", help="file name for log")
    args = parser.parse_args()

    logging.basicConfig(
        stream=sys.stdout,
        level=logging.DEBUG,
        format="%(asctime)s\t%(message)s",
        datefmt="%H:%M:%S",
    )

    logger = logging.getLogger("root")

    if args.log:
        log_file = logging.FileHandler(args.log, mode="w")
        log_file.setFormatter(
            logging.Formatter(fmt="%(asctime)s\t%(message)s", datefmt="%H:%M:%S")
        )
        logger.addHandler(log_file)

    for configuration_file in args.configurations:
        logger.info(os.path.splitext(configuration_file)[0])

        with open(configuration_file) as configuration:
            configurations = json.load(configuration)

        for configuration in configurations:
            logger.info(os.path.splitext(configuration.get("network", ""))[0])

            network = ProteinProteinInteractionNetwork()

            for entry in configuration.get("PTM", {}):
                for gene_name, protein, sites in network.add_proteins_from_excel(
                    entry["file"],
                    entry["label"],
                    entry["time"],
                    protein_accession_col=entry["protein column"],
                    position_col=entry["position column"],
                    replicates=entry["replicate columns"],
                    protein_accession_format=extract.EXTRACT.get(
                        entry["protein format"], lambda x: x
                    ),
                    position_format=extract.EXTRACT.get(
                        entry["position format"], lambda x: x
                    ),
                    sheet_name=entry.get("sheet", 0),
                    header=entry.get("header", 0),
                    num_sites=entry.get("sites", 5),
                    num_replicates=entry.get("replicates", 2),
                    merge_replicates=merge.MERGE.get(
                        entry.get("merge replicates", "mean"), merge.MERGE["mean"]
                    ),
                    convert_measurement=convert.LOG_BASE[entry.get("log base")],
                ):
                    logger.info(
                        "{}\t{}\t{}\t{}\t{}".format(
                            gene_name if gene_name else "",
                            protein,
                            entry["time"],
                            entry["label"],
                            " ".join(["{:.2f}".format(site) for site in sites]),
                        )
                    )

            if configuration.get("PPI"):
                if configuration["PPI"].get("BioGRID"):
                    for (
                        interactor_a,
                        interactor_b,
                    ) in network.add_interactions_from_BioGRID(
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
                            "{}\t{}\tBioGRID\t1.000".format(
                                *sorted([interactor_a, interactor_b])
                            )
                        )

                if configuration["PPI"].get("IntAct"):
                    for (
                        interactor_a,
                        interactor_b,
                        score,
                    ) in network.add_interactions_from_IntAct(
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
                                *sorted([interactor_a, interactor_b]), score
                            )
                        )

                if configuration["PPI"].get("STRING"):
                    for (
                        interactor_a,
                        interactor_b,
                        score,
                    ) in network.add_interactions_from_STRING(
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
                                *sorted([interactor_a, interactor_b]), score
                            )
                        )

            if configuration.get("styles"):
                network.set_post_translational_modification_data_column()

                if configuration.get("Cytoscape", {}).get("type") == "z":
                    network.set_change_data_column(
                        merge_sites=merge.MERGE.get(
                            configuration["Cytoscape"].get("merge sites", "mean"),
                            merge.MERGE["mean"],
                        ),
                        mid_range_thresholds=configuration["Cytoscape"].get(
                            "threshold", (-2.0, 2.0)
                        ),
                        mid_range=network.get_range_z,
                    )
                else:
                    network.set_change_data_column_abs(
                        merge_sites=merge.MERGE.get(
                            configuration["Cytoscape"].get("merge sites", "mean"),
                            merge.MERGE["mean"],
                        ),
                        mid_range_thresholds=configuration["Cytoscape"].get(
                            "threshold", (-1.0, 1.0)
                        ),
                    )

                if configuration.get("network"):
                    export_network(configuration["network"], network)

                styles = CytoscapeStyles(
                    network,
                    bar_chart_range=(
                        configuration.get("Cytoscape", {})
                        .get("bar chart", {})
                        .get("minimum", -3.0),
                        configuration.get("Cytoscape", {})
                        .get("bar chart", {})
                        .get("maximum", 3.0),
                    ),
                )

                export_styles(configuration["styles"], styles)

            if configuration.get("module detection"):
                modules, p_values = network.get_change_enrichment(
                    p=configuration["module detection"].get("p", 0.1),
                    z=configuration["module detection"].get("z", (-2.0, 2.0)),
                    merge_sites=merge.MERGE.get(
                        configuration["module detection"].get("merge sites", "mean"),
                        merge.MERGE["mean"],
                    ),
                    min_size=configuration["module detection"]
                    .get("module size", {})
                    .get("minimum", 3),
                    max_size=configuration["module detection"]
                    .get("module size", {})
                    .get("maximum", 100),
                    merge_sizes=merge.MERGE.get(
                        configuration["module detection"]
                        .get("module size", {})
                        .get("merge", "max"),
                        merge.MERGE["max"],
                    ),
                )
                for time in p_values:
                    for ptm in sorted(p_values[time]):
                        for module in sorted(p_values[time][ptm]):
                            logger.info(
                                "{}\t{}\t{} ({} proteins)\t{:.2E}".format(
                                    time,
                                    ptm,
                                    module + 1,
                                    len(modules[module]),
                                    p_values[time][ptm][module],
                                )
                            )
                            if configuration.get("network"):
                                export_network(
                                    configuration["network"],
                                    network.subgraph(modules[module]),
                                    ".{}".format(module + 1),
                                )


if __name__ == "__main__":
    main()
