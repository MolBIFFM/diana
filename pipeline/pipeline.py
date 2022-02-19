import argparse
import concurrent.futures
import json
import logging
import os
import re
import sys

import networkx as nx

from cytoscape import styles
from databases import biogrid, intact, mint, reactome, string
from interface import combination, conversion, correction, test
from networks import pathway_network, protein_protein_interaction_network
from interface import modularization


def process_configuration(configurations, basename):
    logger = logging.LoggerAdapter(logging.getLogger("root"),
                                   {"configuration": basename})
    for i, configuration in enumerate(configurations, start=1):
        network = nx.Graph()

        for entry in configuration.get("genes", {}):
            if "file" in entry and "accession column" in entry:
                protein_protein_interaction_network.add_genes_from_table(
                    network,
                    file_name=entry["file"],
                    gene_accession_column=entry["accession column"],
                    gene_accession_format=re.compile(
                        entry.get("accession format", "^(.+?)$")),
                    sheet_name=entry.get("sheet", 0),
                    header=entry.get("header", 1) - 1,
                    taxon_identifier=entry.get("taxon identifier", 9606),
                )

            elif "accessions" in entry:
                protein_protein_interaction_network.add_genes_from(
                    network,
                    genes=entry["accessions"],
                    taxon_identifier=entry.get("taxon identifier", 9606),
                )

        for entry in configuration.get("proteins", {}):
            if "file" in entry and "accession column" in entry:
                protein_protein_interaction_network.add_proteins_from_table(
                    network,
                    file_name=entry["file"],
                    protein_accession_column=entry["accession column"],
                    protein_accession_format=re.compile(
                        entry.get("accession format", "^(.+?)$")),
                    time=entry.get("time", 0),
                    modification=entry.get("post-translational modification",
                                           ""),
                    position_column=entry.get("position column", ""),
                    position_format=re.compile(
                        entry.get("position format", "^(.+?)$")),
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
                    measurement_conversion=conversion.LOGARITHM[entry.get(
                        "logarithm")],
                )

            elif "accessions" in entry:
                protein_protein_interaction_network.add_proteins_from(
                    network, proteins=entry["accessions"])

        for entry in configuration.get("networks", []):
            network = nx.algorithms.operators.binary.compose(
                network, nx.readwrite.graphml.read_graphml(entry))

        if "protein-protein interactions" in configuration:
            k = 0
            while any(configuration["protein-protein interactions"].get(
                    database, {}).get("neighbors", 0) > k
                      for database in {"BioGRID", "IntAct", "MINT", "STRING"}):
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

                if "Reactome" in configuration[
                        "protein-protein interactions"] and configuration[
                            "protein-protein interactions"]["Reactome"].get(
                                "neighbors", 0) > k:
                    reactome.add_proteins(
                        network,
                        interaction_context=configuration[
                            "protein-protein interactions"]["Reactome"].get(
                                "interaction context", []),
                        interaction_type=configuration[
                            "protein-protein interactions"]["Reactome"].get(
                                "interaction type", []),
                        taxon_identifier=configuration[
                            "protein-protein interactions"]["Reactome"].get(
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
                        version=configuration["protein-protein interactions"]
                        ["STRING"].get("version", 11.5),
                    )

                k += 1

            if k:
                protein_protein_interaction_network.annotate_proteins(network)
                protein_protein_interaction_network.remove_unannotated_proteins(
                    network)

            if "BioGRID" in configuration["protein-protein interactions"]:
                biogrid.add_protein_protein_interactions(
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

            if "MINT" in configuration["protein-protein interactions"]:
                mint.add_protein_protein_interactions(
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

            if "Reactome" in configuration["protein-protein interactions"]:
                reactome.add_protein_protein_interactions(
                    network,
                    interaction_context=configuration[
                        "protein-protein interactions"]["Reactome"].get(
                            "interaction context", []),
                    interaction_type=configuration[
                        "protein-protein interactions"]["Reactome"].get(
                            "interaction type", []),
                    taxon_identifier=configuration[
                        "protein-protein interactions"]["Reactome"].get(
                            "taxon identifier", 9606),
                )

            if "STRING" in configuration["protein-protein interactions"]:
                string.add_protein_protein_interactions(
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
                    version=configuration["protein-protein interactions"]
                    ["STRING"].get("version", 11.5),
                )

        if "Cytoscape" in configuration:
            if (configuration["Cytoscape"].get("bar chart",
                                               {}).get("type") == "z-score"):
                cytoscape_styles = styles.get_protein_protein_interaction_network_styles(
                    network,
                    bar_chart_range=configuration["Cytoscape"]
                    ["bar chart"].get("range", (-2.0, 2.0)),
                    get_bar_chart_range=protein_protein_interaction_network.
                    get_standard_score_range,
                    site_combination=combination.SITE_COMBINATION.get(
                        configuration["Cytoscape"]["bar chart"].get(
                            "combine sites"),
                        combination.SITE_COMBINATION["absmax"],
                    ),
                    confidence_score_combination=combination.
                    CONFIDENCE_SCORE_COMBINATION.get(
                        configuration["Cytoscape"].get("edge transparency")))

            elif (configuration["Cytoscape"].get(
                    "bar chart", {}).get("type") == "quantile"):
                cytoscape_styles = styles.get_protein_protein_interaction_network_styles(
                    network,
                    bar_chart_range=configuration["Cytoscape"]
                    ["bar chart"].get("range", (0.025, 0.975)),
                    get_bar_chart_range=protein_protein_interaction_network.
                    get_propotion_range,
                    site_combination=combination.SITE_COMBINATION.get(
                        configuration["Cytoscape"]["bar chart"].get(
                            "combine sites"),
                        combination.SITE_COMBINATION["absmax"],
                    ),
                    confidence_score_combination=combination.
                    CONFIDENCE_SCORE_COMBINATION.get(
                        configuration["Cytoscape"].get("edge transparency")))

            else:
                cytoscape_styles = styles.get_protein_protein_interaction_network_styles(
                    network,
                    bar_chart_range=configuration["Cytoscape"]
                    ["bar chart"].get("range", (-1.0, 1.0)),
                    site_combination=combination.SITE_COMBINATION.get(
                        configuration["Cytoscape"]["bar chart"].get(
                            "combine sites"),
                        combination.SITE_COMBINATION["absmax"],
                    ),
                    confidence_score_combination=combination.
                    CONFIDENCE_SCORE_COMBINATION.get(
                        configuration["Cytoscape"].get("edge transparency")))

            protein_protein_interaction_network.set_edge_weights(
                network,
                weight=combination.CONFIDENCE_SCORE_COMBINATION[
                    configuration["Cytoscape"].get("edge transparency")],
                attribute="confidence")

            styles.export(
                cytoscape_styles,
                basename,
                ".{}".format(i) if len(configurations) > 1 else "",
            )

            protein_protein_interaction_network.set_post_translational_modification(
                network)

            if (configuration["Cytoscape"].get("node color",
                                               {}).get("type") == "z-score"):
                protein_protein_interaction_network.set_changes(
                    network,
                    site_combination=combination.SITE_COMBINATION[
                        configuration["Cytoscape"]["node color"].get(
                            "combine sites", "absmax")],
                    changes=configuration["Cytoscape"]["node color"].get(
                        "range", (-2.0, 2.0)),
                    get_range=protein_protein_interaction_network.
                    get_standard_score_range,
                )

            elif (configuration["Cytoscape"].get(
                    "node color", {}).get("type") == "quantile"):
                protein_protein_interaction_network.set_changes(
                    network,
                    site_combination=combination.SITE_COMBINATION[
                        configuration["Cytoscape"]["node color"].get(
                            "combine sites", "absmax")],
                    changes=configuration["Cytoscape"]["node color"].get(
                        "range", (0.025, 0.975)),
                    get_range=protein_protein_interaction_network.
                    get_quantile_range,
                )

            else:
                protein_protein_interaction_network.set_changes(
                    network,
                    site_combination=combination.SITE_COMBINATION[
                        configuration["Cytoscape"]["node color"].get(
                            "combine sites", "absmax")],
                    changes=configuration["Cytoscape"]["node color"].get(
                        "range", (-1.0, 1.0)),
                )

        protein_protein_interaction_network.export(
            network,
            basename,
            ".{}".format(i) if len(configurations) > 1 else "",
        )

        if "post-processing" in configuration:
            if "module detection" in configuration["post-processing"]:
                protein_protein_interaction_network.set_edge_weights(
                    network,
                    weight=combination.CONFIDENCE_SCORE_COMBINATION[
                        configuration["post-processing"]
                        ["module detection"].get("edge weight")])

                for j, module in enumerate(
                        protein_protein_interaction_network.get_modules(
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
                            algorithm=modularization.ALGORITHM.get(
                                configuration["post-processing"]
                                ["module detection"].get("algorithm"),
                                modularization.ALGORITHM["Louvain"],
                            ),
                            resolution=configuration["post-processing"]
                            ["module detection"].get("resolution", 1.0)),
                        start=1,
                ):
                    protein_protein_interaction_network.export(
                        network.subgraph(module),
                        basename,
                        ".{}.{}".format(i, j)
                        if len(configurations) > 1 else ".{}".format(j),
                    )

            if "enrichment analysis" in configuration["post-processing"]:
                protein_protein_interaction_network.set_edge_weights(
                    network,
                    weight=combination.CONFIDENCE_SCORE_COMBINATION[
                        configuration["post-processing"]
                        ["enrichment analysis"].get("edge weight", "number")])

                if (configuration["post-processing"]
                    ["enrichment analysis"].get("type") == "z-score"):
                    modules, p_values = protein_protein_interaction_network.get_change_enriched_modules(
                        network,
                        module_size=configuration["post-processing"]
                        ["enrichment analysis"].get("module size",
                                                    network.number_of_nodes()),
                        p=configuration["post-processing"]
                        ["enrichment analysis"].get("p", 0.05),
                        changes=configuration["post-processing"]
                        ["enrichment analysis"].get("range", (-2.0, 2.0)),
                        get_range=protein_protein_interaction_network.
                        get_standard_score_range,
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
                        algorithm=modularization.ALGORITHM.get(
                            configuration["post-processing"]
                            ["enrichment analysis"].get("algorithm"),
                            modularization.ALGORITHM["Louvain"],
                        ),
                        resolution=configuration["post-processing"]
                        ["enrichment analysis"].get("resolution", 1.0),
                        test=test.TEST.get(
                            configuration["post-processing"]
                            ["enrichment analysis"].get("test"),
                            test.TEST["hypergeometric"],
                        ),
                        correction=correction.CORRECTION.get(
                            configuration["post-processing"]
                            ["enrichment analysis"].get("correction"),
                            correction.CORRECTION["Benjamini-Hochberg"],
                        ))

                elif (configuration["post-processing"]
                      ["enrichment analysis"].get("type") == "quantile"):
                    modules, p_values = protein_protein_interaction_network.get_change_enriched_modules(
                        network,
                        module_size=configuration["post-processing"]
                        ["enrichment analysis"].get("module size",
                                                    network.number_of_nodes()),
                        p=configuration["post-processing"]
                        ["enrichment analysis"].get("p", 0.05),
                        changes=configuration["post-processing"]
                        ["enrichment analysis"].get("range", (0.025, 0.975)),
                        get_range=protein_protein_interaction_network.
                        get_quantile_range,
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
                        algorithm=modularization.ALGORITHM.get(
                            configuration["post-processing"]
                            ["enrichment analysis"].get("algorithm"),
                            modularization.ALGORITHM["Louvain"],
                        ),
                        resolution=configuration["post-processing"]
                        ["enrichment analysis"].get("resolution", 1.0),
                        test=test.TEST.get(
                            configuration["post-processing"]
                            ["enrichment analysis"].get("test"),
                            test.TEST["hypergeometric"],
                        ),
                        correction=correction.CORRECTION.get(
                            configuration["post-processing"]
                            ["enrichment analysis"].get("correction"),
                            correction.CORRECTION["Benjamini-Hochberg"],
                        ))

                else:
                    modules, p_values = protein_protein_interaction_network.get_change_enriched_modules(
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
                        algorithm=modularization.ALGORITHM.get(
                            configuration["post-processing"]
                            ["enrichment analysis"].get("algorithm"),
                            modularization.ALGORITHM["Louvain"],
                        ),
                        resolution=configuration["post-processing"]
                        ["enrichment analysis"].get("resolution", 1.0),
                        test=test.TEST.get(
                            configuration["post-processing"]
                            ["enrichment analysis"].get("test"),
                            test.TEST["hypergeometric"],
                        ),
                        correction=correction.CORRECTION.get(
                            configuration["post-processing"]
                            ["enrichment analysis"].get("correction"),
                            correction.CORRECTION["Benjamini-Hochberg"],
                        ))

                protein_protein_interaction_network.remove_edge_weights(
                    network)

                for time in p_values:
                    for modification in sorted(p_values[time]):
                        for j in sorted(p_values[time][modification]):
                            logger.info("{}\t{}\t{}\t{:.2e}".format(
                                time,
                                modification,
                                j + 1,
                                p_values[time][modification][j],
                            ))
                            protein_protein_interaction_network.export(
                                network.subgraph(modules[j]),
                                basename,
                                ".{}.{}".format(i, j + 1) if
                                len(configurations) > 1 else ".{}".format(j +
                                                                          1),
                            )


def process_configuration_file(configuration_file):
    with open(configuration_file) as configuration:
        process_configuration(
            json.load(configuration),
            os.path.splitext(os.path.basename(configuration_file))[0])


def main():
    logging.basicConfig(stream=sys.stdout,
                        level=logging.INFO,
                        format="%(configuration)s\t%(message)s")

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-c",
        "--configurations",
        help="configuration files",
        nargs="+",
        required=True,
    )
    parser.add_argument(
        "-p",
        "--processes",
        help="maximum number of processes used (default: {})".format(
            os.cpu_count()),
        type=int,
        default=os.cpu_count())

    args = parser.parse_args()

    with concurrent.futures.ProcessPoolExecutor(
            max_workers=args.processes) as executor:
        executor.map(process_configuration_file, args.configurations)


if __name__ == "__main__":
    main()
