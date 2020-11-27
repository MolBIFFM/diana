import argparse
import logging

import yaml

from pipeline import protein_protein_interaction_network
from pipeline import cytoscape_style


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-c",
                        "--configuration",
                        help="YAML configuration file")
    parser.add_argument("-n",
                        "--network",
                        help="protein-protein interaction network")
    parser.add_argument("-s", "--styles", help="Cytoscape style specification")
    parser.add_argument("-l", "--log", help="log")
    args = parser.parse_args()

    with open(args.configuration) as configuration_file:
        configuration = yaml.load(configuration_file, Loader=yaml.Loader)

    if args.log:
        logging.basicConfig(filename=args.log,
                            filemode="w",
                            level=logging.DEBUG,
                            format="%(asctime)s\t%(levelname)s\t%(message)s",
                            datefmt="%H:%M:%S")

    ppi_network = protein_protein_interaction_network.ProteinProteinInteractionNetwork(
    )

    for entry in configuration.get("PTM", {}):
        for protein in ppi_network.add_proteins_from_excel(
                entry["file"],
                entry["label"],
                entry["time"],
                num_replicates=entry.get("replicates", 2),
                num_sites=entry.get("sites", 5)):
            if args.log:
                logging.info("{} {} {}".format(entry["label"], entry["time"],
                                               protein))

    if configuration.get("PPI"):
        if configuration["PPI"].get("BioGRID"):
            for interactor_a, interactor_b in ppi_network.add_interactions_from_BioGRID(
                    experimental_system=configuration["PPI"]
                ["BioGRID"].get("experimental system", [
                    "Affinity Capture-Luminescence", "Affinity Capture-MS",
                    "Affinity Capture-RNA", "Affinity Capture-Western",
                    "Biochemical Activity", "Co-crystal Structure",
                    "Co-purification", "FRET", "PCA", "Two-hybrid"
                ])):
                if args.log:
                    logging.info("BioGRID: {}, {} 1.0".format(
                        interactor_a, interactor_b))

        if configuration["PPI"].get("IntAct"):
            for interactor_a, interactor_b, score in ppi_network.add_interactions_from_IntAct(
                    interaction_detection_methods=configuration["PPI"]
                ["IntAct"].get("interaction detection methods", [
                    "affinity chromatography technology", "two hybrid",
                    "biochemical", "pull down", "enzymatic study", "bio id",
                    "x-ray crystallography",
                    "fluorescent resonance energy transfer",
                    "protein complementation assay"
                ]),
                    interaction_types=configuration["PPI"]["IntAct"].get(
                        "interaction_types", [
                            "physical association", "direct interaction",
                            "association"
                        ]),
                    mi_score=configuration["PPI"]["IntAct"].get(
                        "MI score", 0.27)):
                if args.log:
                    logging.info("IntAct: {}, {} {}".format(
                        interactor_a, interactor_b, score))

        if configuration["PPI"].get("STRING"):
            for interactor_a, interactor_b, score in ppi_network.add_interactions_from_STRING(
                    neighborhood=configuration["PPI"]["STRING"].get(
                        "neighborhood", 0.0),
                    neighborhood_transferred=configuration["PPI"]
                ["STRING"].get("neighborhood transferred", 0.0),
                    fusion=configuration["PPI"]["STRING"].get("fusion", 0.0),
                    cooccurence=configuration["PPI"]["STRING"].get(
                        "cooccurence", 0.0),
                    homology=configuration["PPI"]["STRING"].get(
                        "homology", 0.0),
                    coexpression=configuration["PPI"]["STRING"].get(
                        "coexpression", 0.0),
                    coexpression_transferred=configuration["PPI"]
                ["STRING"].get("coexpression transferred", 0.0),
                    experiments=configuration["PPI"]["STRING"].get(
                        "experiments", 0.7),
                    experiments_transferred=configuration["PPI"]["STRING"].get(
                        "experiments transferred", 0.0),
                    database=configuration["PPI"]["STRING"].get(
                        "database", 0.0),
                    database_transferred=configuration["PPI"]["STRING"].get(
                        "database transferred", 0.0),
                    textmining=configuration["PPI"]["STRING"].get(
                        "textmining", 0.0),
                    textmining_transferred=configuration["PPI"]["STRING"].get(
                        "textmining transferred", 0.0),
                    combined_score=configuration["PPI"]["STRING"].get(
                        "combined score", 0.7)):
                if args.log:
                    logging.info("STRING {}, {} {}".format(
                        interactor_a, interactor_b, score))

    ppi_network.remove_isolates()

    if args.styles:
        style = cytoscape_style.CytoscapeStyle(
            ppi_network,
            bar_chart_range=(configuration.get("Cytoscape", {}).get(
                "bar chart", {}).get("minimum", -3.0),
                             configuration.get("Cytoscape",
                                               {}).get("bar chart",
                                                       {}).get("maximum",
                                                               3.0)))

        ppi_network.set_ptm_data_column()
        ppi_network.set_trend_data_column(mid_range=(
            configuration.get("Cytoscape", {}).get("mid range", {}).get(
                "minimum", -1.0),
            configuration.get("Cytoscape", {}).get("mid range", {}).get(
                "maximum", 1.0),
        ))

        if args.styles.endswith(".xml"):
            style.export_as_xml(args.styles)

    if args.network:
        if args.network.endswith(".graphml") or args.network.endswith(".xml"):
            ppi_network.export_as_graphml(args.network)
        elif args.network.endswith(".cyjs") or args.network.endswith(".json"):
            ppi_network.export_as_cyjs(args.network)


if __name__ == "__main__":
    main()
