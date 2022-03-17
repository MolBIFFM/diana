"""pipeline"""
import argparse
import concurrent.futures
import json
import logging
import os
import re
import sys
from typing import Optional

import networkx as nx

from cytoscape import styles
from interface import combination, conversion, correction, modularization, test
from networks import (ontology_network, pathway_network,
                      protein_protein_interaction_network)


def process_configuration(configurations: list[dict],
                          logger: logging.Logger) -> None:
    """
    Executes workflows specified in configurations sequentially.

    Args:
        configurations: The specification of workflows.
        logger: A configuration-specific logger.
    """
    for i, configuration in enumerate(configurations, start=1):
        process_workflow(configuration, logger,
                         i if len(configurations) > 1 else None)


def process_workflow(configuration: dict,
                     logger: logging.Logger,
                     index: Optional[int] = None) -> None:
    """
        Executes a workflow specified in configuration.

        Args:
            configurations: The specification of a workflow.
            logger: A logger.
    """
    network = protein_protein_interaction_network.get_protein_protein_interaction_network(
    )

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
                taxonomy_identifier=entry.get("taxonomy identifier", 9606),
            )

        elif "accessions" in entry:
            protein_protein_interaction_network.add_genes_from(
                network,
                genes=entry["accessions"],
                taxonomy_identifier=entry.get("taxonomy identifier", 9606),
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
                modification=entry.get("post-translational modification", ""),
                position_column=entry.get("position column", ""),
                position_format=re.compile(
                    entry.get("position format", "^(.+?)$")),
                replicates=entry.get("replicate columns", []),
                sheet_name=entry.get("sheet", 0),
                header=entry.get("header", 1) - 1,
                num_sites=entry.get("sites", 5),
                num_replicates=entry.get("replicates", 1),
                replicate_combination=combination.REPLICATE_COMBINATION[
                    entry.get("combine replicates", "mean")],
                measurement_conversion=conversion.LOGARITHM[entry.get(
                    "logarithm")])

        elif "accessions" in entry:
            protein_protein_interaction_network.add_proteins_from(
                network, proteins=entry["accessions"])

    for entry in configuration.get("networks", []):
        network = nx.compose(network, nx.readwrite.graphml.read_graphml(entry))

    if "protein-protein interactions" in configuration:
        neighbors = 0
        while any(configuration["protein-protein interactions"].get(
                database, {}).get("neighbors", 0) > neighbors for database in
                  {"BioGRID", "IntAct", "MINT", "Reactome", "STRING"}):
            if "BioGRID" in configuration[
                    "protein-protein interactions"] and configuration[
                        "protein-protein interactions"]["BioGRID"].get(
                            "neighbors", 0) > neighbors:
                protein_protein_interaction_network.add_proteins_from_biogrid(
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
                    taxonomy_identifier=configuration[
                        "protein-protein interactions"]["BioGRID"].get(
                            "taxonomy identifier", 9606),
                    version=tuple(
                        int(v)
                        for v in configuration["protein-protein interactions"]
                        ["BioGRID"].get("version").split(".")) if "version"
                    in configuration["protein-protein interactions"]["BioGRID"]
                    else None)

            if "IntAct" in configuration[
                    "protein-protein interactions"] and configuration[
                        "protein-protein interactions"]["IntAct"].get(
                            "neighbors", 0) > neighbors:
                protein_protein_interaction_network.add_proteins_from_intact(
                    network,
                    interaction_detection_methods=configuration[
                        "protein-protein interactions"]["IntAct"].get(
                            "interaction detection methods", []),
                    interaction_types=configuration[
                        "protein-protein interactions"]["IntAct"].get(
                            "interaction types", []),
                    mi_score=configuration["protein-protein interactions"]
                    ["IntAct"].get("MI score", 0.0),
                    taxonomy_identifier=configuration[
                        "protein-protein interactions"]["IntAct"].get(
                            "taxonomy identifier", 9606),
                )

            if "MINT" in configuration[
                    "protein-protein interactions"] and configuration[
                        "protein-protein interactions"]["MINT"].get(
                            "neighbors", 0) > neighbors:
                protein_protein_interaction_network.add_proteins_from_mint(
                    network,
                    interaction_detection_methods=configuration[
                        "protein-protein interactions"]["MINT"].get(
                            "interaction detection methods", []),
                    interaction_types=configuration[
                        "protein-protein interactions"]["MINT"].get(
                            "interaction types", []),
                    mi_score=configuration["protein-protein interactions"]
                    ["MINT"].get("MI score", 0.0),
                    taxonomy_identifier=configuration[
                        "protein-protein interactions"]["MINT"].get(
                            "taxonomy identifier", 9606),
                )

            if "Reactome" in configuration[
                    "protein-protein interactions"] and configuration[
                        "protein-protein interactions"]["Reactome"].get(
                            "neighbors", 0) > neighbors:
                protein_protein_interaction_network.add_proteins_from_reactome(
                    network,
                    interaction_context=configuration[
                        "protein-protein interactions"]["Reactome"].get(
                            "interaction context", []),
                    interaction_type=configuration[
                        "protein-protein interactions"]["Reactome"].get(
                            "interaction type", []),
                    taxonomy_identifier=configuration[
                        "protein-protein interactions"]["Reactome"].get(
                            "taxonomy identifier", 9606),
                )

            if "STRING" in configuration[
                    "protein-protein interactions"] and configuration[
                        "protein-protein interactions"]["STRING"].get(
                            "neighbors", 0) > neighbors:
                protein_protein_interaction_network.add_proteins_from_string(
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
                    combined_score=configuration["protein-protein interactions"]
                    ["STRING"].get("combined score", 0.0),
                    physical=configuration["protein-protein interactions"]
                    ["STRING"].get("physical", False),
                    taxonomy_identifier=configuration[
                        "protein-protein interactions"]["STRING"].get(
                            "taxonomy identifier", 9606),
                    version=configuration["protein-protein interactions"]
                    ["STRING"].get("version", 11.5),
                )

            neighbors += 1

        if neighbors:
            for taxonomy_identifier in (
                    configuration["protein-protein interactions"][database].get(
                        "taxonomy_identifier", 9606) for database in
                {"BioGRID", "IntAct", "MINT", "Reactome", "STRING"}):
                protein_protein_interaction_network.annotate_proteins(
                    network, taxonomy_identifier)

            protein_protein_interaction_network.remove_unannotated_proteins(
                network)

        if "BioGRID" in configuration["protein-protein interactions"]:
            protein_protein_interaction_network.add_protein_protein_interactions_from_biogrid(
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
                taxonomy_identifier=configuration[
                    "protein-protein interactions"]["BioGRID"].get(
                        "taxonomy identifier", 9606),
                version=tuple(
                    int(v)
                    for v in configuration["protein-protein interactions"]
                    ["BioGRID"].get("version").split(".")) if "version"
                in configuration["protein-protein interactions"]["BioGRID"] else
                None)

        if "IntAct" in configuration["protein-protein interactions"]:
            protein_protein_interaction_network.add_protein_protein_interactions_from_intact(
                network,
                interaction_detection_methods=configuration[
                    "protein-protein interactions"]["IntAct"].get(
                        "interaction detection methods", []),
                interaction_types=configuration["protein-protein interactions"]
                ["IntAct"].get("interaction types", []),
                mi_score=configuration["protein-protein interactions"]
                ["IntAct"].get("MI score", 0.0),
                taxonomy_identifier=configuration[
                    "protein-protein interactions"]["IntAct"].get(
                        "taxonomy identifier", 9606),
            )

        if "MINT" in configuration["protein-protein interactions"]:
            protein_protein_interaction_network.add_protein_protein_interactions_from_mint(
                network,
                interaction_detection_methods=configuration[
                    "protein-protein interactions"]["MINT"].get(
                        "interaction detection methods", []),
                interaction_types=configuration["protein-protein interactions"]
                ["MINT"].get("interaction types", []),
                mi_score=configuration["protein-protein interactions"]
                ["MINT"].get("MI score", 0.0),
                taxonomy_identifier=configuration[
                    "protein-protein interactions"]["MINT"].get(
                        "taxonomy identifier", 9606),
            )

        if "Reactome" in configuration["protein-protein interactions"]:
            protein_protein_interaction_network.add_protein_protein_interactions_from_reactome(
                network,
                interaction_context=configuration[
                    "protein-protein interactions"]["Reactome"].get(
                        "interaction context", []),
                interaction_type=configuration["protein-protein interactions"]
                ["Reactome"].get("interaction type", []),
                taxonomy_identifier=configuration[
                    "protein-protein interactions"]["Reactome"].get(
                        "taxonomy identifier", 9606),
            )

        if "STRING" in configuration["protein-protein interactions"]:
            protein_protein_interaction_network.add_protein_protein_interactions_from_string(
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
                combined_score=configuration["protein-protein interactions"]
                ["STRING"].get("combined score", 0.0),
                physical=configuration["protein-protein interactions"]
                ["STRING"].get("physical", False),
                taxonomy_identifier=configuration[
                    "protein-protein interactions"]["STRING"].get(
                        "taxonomy identifier", 9606),
                version=configuration["protein-protein interactions"]
                ["STRING"].get("version", 11.5),
            )

    if "Cytoscape" in configuration:
        if (configuration["Cytoscape"].get(
                "bar chart", {}).get("conversion") == "standard score"):
            cytoscape_styles = styles.get_protein_protein_interaction_network_styles(
                network,
                bar_chart_range=configuration["Cytoscape"]["bar chart"].get(
                    "range", (-2.0, 2.0)),
                convert_change=protein_protein_interaction_network.
                convert_standard_score_to_change,
                site_combination=combination.SITE_COMBINATION[
                    configuration["Cytoscape"]["bar chart"].get(
                        "site combination", "absmax")],
                confidence_score_combination=combination.
                CONFIDENCE_SCORE_COMBINATION[configuration["Cytoscape"].get(
                    "edge transparency")])

        elif (configuration["Cytoscape"].get(
                "bar chart", {}).get("conversion") == "quantile"):
            cytoscape_styles = styles.get_protein_protein_interaction_network_styles(
                network,
                bar_chart_range=configuration["Cytoscape"]["bar chart"].get(
                    "range", (0.025, 0.975)),
                convert_change=protein_protein_interaction_network.
                convert_quantile_to_change,
                site_combination=combination.SITE_COMBINATION[
                    configuration["Cytoscape"]["bar chart"].get(
                        "site combination", "absmax")],
                confidence_score_combination=combination.
                CONFIDENCE_SCORE_COMBINATION[configuration["Cytoscape"].get(
                    "edge transparency")])

        else:
            cytoscape_styles = styles.get_protein_protein_interaction_network_styles(
                network,
                bar_chart_range=configuration["Cytoscape"]["bar chart"].get(
                    "range", (-1.0, 1.0)),
                site_combination=combination.SITE_COMBINATION[
                    configuration["Cytoscape"]["bar chart"].get(
                        "site combination", "absmax")],
                confidence_score_combination=combination.
                CONFIDENCE_SCORE_COMBINATION[configuration["Cytoscape"].get(
                    "edge transparency")])

        protein_protein_interaction_network.set_edge_weights(
            network,
            weight=combination.CONFIDENCE_SCORE_COMBINATION[
                configuration["Cytoscape"].get("edge transparency")],
            attribute="score")

        styles.export(
            cytoscape_styles,
            logger.name,
            ".{}".format(index) if index else "",
        )

        protein_protein_interaction_network.set_post_translational_modification(
            network)

        if (configuration["Cytoscape"].get(
                "node color", {}).get("conversion") == "standard score"):
            protein_protein_interaction_network.set_changes(
                network,
                site_combination=combination.SITE_COMBINATION[
                    configuration["Cytoscape"]["node color"].get(
                        "site combination", "absmax")],
                changes=configuration["Cytoscape"]["node color"].get(
                    "change", (-2.0, 2.0)),
                convert_change=protein_protein_interaction_network.
                convert_standard_score_to_change,
            )

        elif (configuration["Cytoscape"].get(
                "node color", {}).get("conversion") == "quantile"):
            protein_protein_interaction_network.set_changes(
                network,
                site_combination=combination.SITE_COMBINATION[
                    configuration["Cytoscape"]["node color"].get(
                        "site combination", "absmax")],
                changes=configuration["Cytoscape"]["node color"].get(
                    "change", (0.025, 0.975)),
                convert_change=protein_protein_interaction_network.
                convert_quantile_to_change,
            )

        else:
            protein_protein_interaction_network.set_changes(
                network,
                site_combination=combination.SITE_COMBINATION[
                    configuration["Cytoscape"]["node color"].get(
                        "site combination", "absmax")],
                changes=configuration["Cytoscape"]["node color"].get(
                    "change", (-1.0, 1.0)),
            )

    protein_protein_interaction_network.export(
        network,
        logger.name,
        ".{}".format(index) if index else "",
    )

    if "Gene Ontology enrichment" in configuration:
        enrichment = protein_protein_interaction_network.get_gene_ontology_enrichment(
            [network],
            enrichment_test=test.ENRICHMENT_TEST[
                configuration["module detection"]
                ["Gene Ontology enrichment"].get("test", "hypergeometric")],
            multiple_testing_correction=correction.CORRECTION[
                configuration["module detection"].get("correction",
                                                      "Benjamini-Hochberg")],
            taxonomy_identifier=configuration["module detection"]
            ["Gene Ontology enrichment"].get("taxonomy identifier", 9606),
            namespaces=configuration["module detection"]
            ["Gene Ontology enrichment"].get("namespaces", [
                "cellular_component", "molecular_function", "biological_process"
            ]))

        for (term, name), p in sorted(enrichment[network].items(),
                                      key=lambda item: item[1]):
            if p <= configuration["Gene Ontology enrichment"].get("p", 1.0):
                logger.info("{}\t{:.2e}\t{}".format(
                    term,
                    p,
                    name,
                ))

    if "module detection" in configuration:
        protein_protein_interaction_network.set_edge_weights(
            network,
            weight=combination.CONFIDENCE_SCORE_COMBINATION[
                configuration["module detection"].get("edge weight")])

        modules = protein_protein_interaction_network.get_modules(
            network,
            module_size=configuration["module detection"].get(
                "module size", network.number_of_nodes()),
            module_size_combination=combination.MODULE_SIZE_COMBINATION[
                configuration["module detection"].get("module size combination",
                                                      "mean")],
            algorithm=modularization.ALGORITHM[
                configuration["module detection"].get("algorithm", "Louvain")],
            resolution=configuration["module detection"].get("resolution", 1.0))

        protein_protein_interaction_network.remove_edge_weights(network)

        if "Gene Ontology enrichment" in configuration["module detection"]:
            gene_ontology_enrichment = {}
            if "annotation" in configuration["module detection"][
                    "Gene Ontology enrichment"]:
                gene_ontology_enrichment[
                    "annotation"] = protein_protein_interaction_network.get_gene_ontology_enrichment(
                        modules,
                        enrichment_test=test.ENRICHMENT_TEST[
                            configuration["module detection"]
                            ["Gene Ontology enrichment"]["annotation"].get(
                                "test", "hypergeometric")],
                        multiple_testing_correction=correction.CORRECTION[
                            configuration["module detection"]
                            ["Gene Ontology enrichment"]["annotation"].get(
                                "correction", "Benjamini-Hochberg")],
                        taxonomy_identifier=configuration["module detection"]
                        ["Gene Ontology enrichment"]["annotation"].get(
                            "taxonomy identifier", 9606),
                        namespaces=configuration["module detection"]
                        ["Gene Ontology enrichment"]["annotation"].get(
                            "namespaces", [
                                "cellular_component", "molecular_function",
                                "biological_process"
                            ]))
            if "network" in configuration["module detection"][
                    "Gene Ontology enrichment"]:
                gene_ontology_enrichment[
                    "network"] = protein_protein_interaction_network.get_gene_ontology_enrichment(
                        modules,
                        enrichment_test=test.ENRICHMENT_TEST[
                            configuration["module detection"]
                            ["Gene Ontology enrichment"]["network"].get(
                                "test", "hypergeometric")],
                        multiple_testing_correction=correction.CORRECTION[
                            configuration["module detection"]
                            ["Gene Ontology enrichment"]["network"].get(
                                "correction", "Benjamini-Hochberg")],
                        taxonomy_identifier=configuration["module detection"]
                        ["Gene Ontology enrichment"]["network"].get(
                            "taxonomy identifier", 9606),
                        namespaces=configuration["module detection"]
                        ["Gene Ontology enrichment"]["network"].get(
                            "namespaces", [
                                "cellular_component", "molecular_function",
                                "biological_process"
                            ]),
                        annotation_as_reference=False)

        if "change enrichment" in configuration["module detection"]:
            change_enrichment = {}
            if "proteins" in configuration["module detection"][
                    "change enrichment"]:
                if (configuration["module detection"]["change enrichment"]
                    ["proteins"].get("conversion") == "standard score"):
                    change_enrichment[
                        "proteins"] = protein_protein_interaction_network.get_change_enrichment(
                            network,
                            modules,
                            changes=configuration["module detection"]
                            ["change enrichment"]["proteins"].get(
                                "change", (-2.0, 2.0)),
                            convert_change=protein_protein_interaction_network.
                            convert_standard_score_to_change,
                            site_combination=combination.SITE_COMBINATION[
                                configuration["module detection"]
                                ["change enrichment"]["proteins"].get(
                                    "site combination", "absmax")],
                            enrichment_test=test.ENRICHMENT_TEST[
                                configuration["module detection"]
                                ["change enrichment"]["proteins"].get(
                                    "test", "hypergeometric")],
                            multiple_testing_correction=correction.CORRECTION[
                                configuration["module detection"]
                                ["change enrichment"]["proteins"].get(
                                    "correction", "Benjamini-Hochberg")])

                elif (configuration["module detection"]["change enrichment"]
                      ["proteins"].get("conversion") == "quantile"):
                    change_enrichment[
                        "proteins"] = protein_protein_interaction_network.get_change_enrichment(
                            network,
                            modules,
                            changes=configuration["module detection"]
                            ["change enrichment"]["proteins"].get(
                                "change", (0.025, 0.975)),
                            convert_change=protein_protein_interaction_network.
                            convert_quantile_to_change,
                            site_combination=combination.SITE_COMBINATION[
                                configuration["module detection"]
                                ["change enrichment"]["proteins"].get(
                                    "site combination", "absmax")],
                            enrichment_test=test.ENRICHMENT_TEST[
                                configuration["module detection"]
                                ["change enrichment"]["proteins"].get(
                                    "test", "hypergeometric")],
                            multiple_testing_correction=correction.CORRECTION[
                                configuration["module detection"]
                                ["change enrichment"]["proteins"].get(
                                    "correction", "Benjamini-Hochberg")],
                        )

                else:
                    change_enrichment[
                        "proteins"] = protein_protein_interaction_network.get_change_enrichment(
                            network,
                            modules,
                            changes=configuration["module detection"]
                            ["change enrichment"]["proteins"].get(
                                "change", (-1.0, 1.0)),
                            site_combination=combination.SITE_COMBINATION[
                                configuration["module detection"]
                                ["change enrichment"]["proteins"].get(
                                    "site combination", "absmax")],
                            enrichment_test=test.ENRICHMENT_TEST[
                                configuration["module detection"]
                                ["change enrichment"]["proteins"].get(
                                    "test", "hypergeometric")],
                            multiple_testing_correction=correction.CORRECTION[
                                configuration["module detection"]
                                ["change enrichment"]["proteins"].get(
                                    "correction", "Benjamini-Hochberg")],
                        )

            if "sites" in configuration["module detection"][
                    "change enrichment"]:
                if (configuration["module detection"]["change enrichment"]
                    ["sites"].get("conversion") == "standard score"):
                    change_enrichment[
                        "sites"] = protein_protein_interaction_network.get_change_enrichment(
                            network,
                            modules,
                            changes=configuration["module detection"]
                            ["change enrichment"]["proteins sites"].get(
                                "change", (-2.0, 2.0)),
                            convert_change=protein_protein_interaction_network.
                            convert_standard_score_to_change,
                            enrichment_test=test.ENRICHMENT_TEST[
                                configuration["module detection"]
                                ["change enrichment"]["sites"].get(
                                    "test", "hypergeometric")],
                            multiple_testing_correction=correction.CORRECTION[
                                configuration["module detection"]
                                ["change enrichment"]["sites"].get(
                                    "correction", "Benjamini-Hochberg")])

                elif (configuration["module detection"]["change enrichment"]
                      ["sites"].get("conversion") == "quantile"):
                    change_enrichment[
                        "sites"] = protein_protein_interaction_network.get_change_enrichment(
                            network,
                            modules,
                            changes=configuration["module detection"]
                            ["change enrichment"]["sites"].get(
                                "change", (0.025, 0.975)),
                            convert_change=protein_protein_interaction_network.
                            convert_quantile_to_change,
                            enrichment_test=test.ENRICHMENT_TEST[
                                configuration["module detection"]
                                ["change enrichment"]["sites"].get(
                                    "test", "hypergeometric")],
                            multiple_testing_correction=correction.CORRECTION[
                                configuration["module detection"]
                                ["change enrichment"]["sites"].get(
                                    "correction", "Benjamini-Hochberg")],
                        )

                else:
                    change_enrichment[
                        "sites"] = protein_protein_interaction_network.get_change_enrichment(
                            network,
                            modules,
                            changes=configuration["module detection"]
                            ["change enrichment"]["sites"].get(
                                "change", (-1.0, 1.0)),
                            enrichment_test=test.ENRICHMENT_TEST[
                                configuration["module detection"]
                                ["change enrichment"]["sites"].get(
                                    "test", "hypergeometric")],
                            multiple_testing_correction=correction.CORRECTION[
                                configuration["module detection"]
                                ["change enrichment"]["sites"].get(
                                    "correction", "Benjamini-Hochberg")],
                        )
        if "change location" in configuration["module detection"]:
            change_location = {}
            if "proteins" in configuration["module detection"][
                    "change location"]:
                change_location[
                    "proteins"] = protein_protein_interaction_network.get_change_location(
                        network,
                        modules,
                        combination.SITE_COMBINATION[
                            configuration["module detection"]["change location"]
                            ["proteins"].get("site combination", "absmax")],
                        location_test=test.LOCATION_TEST[
                            configuration["module detection"]["change location"]
                            ["proteins"].get("test", "Wilcoxon")],
                        multiple_testing_correction=correction.CORRECTION[
                            configuration["module detection"]["change location"]
                            ["proteins"].get("correction",
                                             "Benjamini-Hochberg")])
            if "sites" in configuration["module detection"]["change location"]:
                change_location[
                    "sites"] = protein_protein_interaction_network.get_change_location(
                        network,
                        modules,
                        location_test=test.LOCATION_TEST[
                            configuration["module detection"]["change location"]
                            ["sites"].get("test", "Wilcoxon")],
                        multiple_testing_correction=correction.CORRECTION[
                            configuration["module detection"]["change location"]
                            ["sites"].get("correction", "Benjamini-Hochberg")])

        for k, module in enumerate(
                sorted(modules,
                       key=lambda module: module.number_of_nodes(),
                       reverse=True),
                start=1):
            export = False
            if "Gene Ontology enrichment" in configuration["module detection"]:
                if "annotation" in configuration["module detection"][
                        "Gene Ontology enrichment"]:
                    for (term, name), p in sorted(
                            gene_ontology_enrichment["annotation"]
                        [module].items(),
                            key=lambda item: item[1]):
                        if p <= configuration["module detection"][
                                "Gene Ontology enrichment"]["annotation"].get(
                                    "p", 1.0):
                            export = True
                            logger.info("{} annotation\t{}\t{:.2e}\t{}".format(
                                k, term, p, name))

                if "network" in configuration["module detection"][
                        "Gene Ontology enrichment"]:
                    for (term, name), p in sorted(
                            gene_ontology_enrichment["network"][module].items(),
                            key=lambda item: item[1]):
                        if p <= configuration["module detection"][
                                "Gene Ontology enrichment"]["network"].get(
                                    "p", 1.0):
                            export = True
                            logger.info("{} network\t{}\t{:.2e}\t{}".format(
                                k, term, p, name))

            if "change enrichment" in configuration["module detection"]:
                if "proteins" in configuration["module detection"][
                        "change enrichment"]:
                    for time in change_enrichment["proteins"][module]:
                        for modification, p in sorted(
                                change_enrichment["proteins"][module]
                            [time].items(),
                                key=lambda item: item[1]):
                            if p <= configuration["module detection"][
                                    "change enrichment"]["proteins"].get(
                                        "p", 1.0):
                                export = True
                                logger.info(
                                    "{}\tprotein change enrichment\t{} {}\t{:.2e}"
                                    .format(
                                        k,
                                        time,
                                        modification,
                                        p,
                                    ))

                if "sites" in configuration["module detection"][
                        "change enrichment"]:
                    for time in change_enrichment["sites"][module]:
                        for modification, p in sorted(change_enrichment["sites"]
                                                      [module][time].items(),
                                                      key=lambda item: item[1]):
                            if p <= configuration["module detection"][
                                    "change enrichment"]["sites"].get("p", 1.0):
                                export = True
                                logger.info(
                                    "{}\tsite change enrichment\t{} {}\t{:.2e}".
                                    format(
                                        k,
                                        time,
                                        modification,
                                        p,
                                    ))

            if "change location" in configuration["module detection"]:
                if "proteins" in configuration["module detection"][
                        "change location"]:
                    for time in change_location["proteins"][module]:
                        for modification, p in sorted(
                                change_location["proteins"][module]
                            [time].items(),
                                key=lambda item: item[1]):
                            if p <= configuration["module detection"][
                                    "change location"]["proteins"].get(
                                        "p", 1.0):
                                export = True
                                logger.info(
                                    "{}\tprotein change location\t{} {}\t{:.2e}"
                                    .format(
                                        k,
                                        time,
                                        modification,
                                        p,
                                    ))
                if "sites" in configuration["module detection"][
                        "change location"]:
                    for time in change_location["sites"][module]:
                        for modification, p in sorted(
                                change_location["sites"][module][time].items(),
                                key=lambda item: item[1]):
                            if p <= configuration["module detection"][
                                    "change location"]["sites"].get("p", 1.0):
                                export = True
                                logger.info(
                                    "{}\tsite change location\t{} {}\t{:.2e}".
                                    format(
                                        k,
                                        time,
                                        modification,
                                        p,
                                    ))

            if export:
                protein_protein_interaction_network.export(
                    module,
                    logger.name,
                    ".{}.{}".format(index, k) if index else ".{}".format(k),
                )
    if "Reactome network" in configuration:
        if "union" in configuration["Reactome network"]:
            proteins = set()
            for subset in configuration["Reactome network"]["union"]:
                if subset.get(
                        "time"
                ) in protein_protein_interaction_network.get_times(
                        network
                ) and subset.get(
                        "post-translational modification"
                ) in protein_protein_interaction_network.get_post_translational_modifications(
                        network, subset["time"]):
                    if subset.get("conversion") == "standard score":
                        change_range = (
                            protein_protein_interaction_network.
                            convert_standard_score_to_change(
                                network,
                                subset.get("change",
                                           (-2.0, 2.0))[0], subset["time"],
                                subset["post-translational modification"],
                                combination.SITE_COMBINATION[subset.get(
                                    "site combination", "absmax")]),
                            protein_protein_interaction_network.
                            convert_standard_score_to_change(
                                network,
                                subset.get("change",
                                           (-2.0, 2.0))[1], subset["time"],
                                subset["post-translational modification"],
                                combination.SITE_COMBINATION[subset.get(
                                    "site combination", "absmax")]))
                        proteins.update(
                            protein_protein_interaction_network.get_proteins(
                                network, subset["time"],
                                subset["post-translational modification"],
                                combination.SITE_COMBINATION[subset.get(
                                    "site combination", "absmax")],
                                lambda change: change <= change_range[
                                    0] or change >= change_range[1]))

                    elif subset.get("conversion") == "quantile":
                        change_range = (
                            protein_protein_interaction_network.
                            convert_quantile_to_change(
                                network,
                                subset.get("change",
                                           (0.025, 0.975))[0], subset["time"],
                                subset["post-translational modification"],
                                combination.SITE_COMBINATION[subset.get(
                                    "site combination", "absmax")]),
                            protein_protein_interaction_network.
                            convert_quantile_to_change(
                                network,
                                subset.get("change",
                                           (0.025, 0.975))[1], subset["time"],
                                subset["post-translational modification"],
                                combination.SITE_COMBINATION[subset.get(
                                    "site combination", "absmax")]))
                        proteins.update(
                            protein_protein_interaction_network.get_proteins(
                                network, subset["time"],
                                subset["post-translational modification"],
                                combination.SITE_COMBINATION[subset.get(
                                    "site combination", "absmax")],
                                lambda change: change <= change_range[
                                    0] or change >= change_range[1]))
                    else:
                        change_range = (subset.get("change", (-1.0, 1.0))[0],
                                        subset.get("change", (-1.0, 1.0))[1])
                        proteins.update(
                            protein_protein_interaction_network.get_proteins(
                                network, subset["time"],
                                subset["post-translational modification"],
                                combination.SITE_COMBINATION[subset.get(
                                    "site combination", "absmax")],
                                lambda change: change <= change_range[
                                    0] or change >= change_range[1]))

            reactome_network = pathway_network.get_pathway_network(
                nx.induced_subgraph(network, proteins),
                enrichment_test=test.ENRICHMENT_TEST[
                    configuration["Reactome network"].get(
                        "test", "hypergeometric")],
                multiple_testing_correction=correction.CORRECTION[
                    configuration["Reactome network"].get(
                        "correction", "Benjamini-Hochberg")],
                taxonomy_identifier=configuration["Reactome network"].get(
                    "taxonomy identifier", 9606))

        elif "intersection" in configuration["Reactome network"]:
            proteins = set(network.nodes())
            for subset in configuration["Reactome network"]["intersection"]:
                if subset.get(
                        "time"
                ) in protein_protein_interaction_network.get_times(
                        network
                ) and subset.get(
                        "post-translational modification"
                ) in protein_protein_interaction_network.get_post_translational_modifications(
                        network, subset["time"]):
                    if subset.get("conversion") == "standard score":
                        change_range = (
                            protein_protein_interaction_network.
                            convert_standard_score_to_change(
                                network,
                                subset.get("change",
                                           (-2.0, 2.0))[0], subset["time"],
                                subset["post-translational modification"],
                                combination.SITE_COMBINATION[subset.get(
                                    "site combination", "absmax")]),
                            protein_protein_interaction_network.
                            convert_standard_score_to_change(
                                network,
                                subset.get("change",
                                           (-2.0, 2.0))[1], subset["time"],
                                subset["post-translational modification"],
                                combination.SITE_COMBINATION[subset.get(
                                    "site combination", "absmax")]))
                        proteins.intersection_update(
                            protein_protein_interaction_network.get_proteins(
                                network, subset["time"],
                                subset["post-translational modification"],
                                combination.SITE_COMBINATION[subset.get(
                                    "site combination", "absmax")],
                                lambda change: change <= change_range[
                                    0] or change >= change_range[1]))

                    elif subset.get("conversion") == "quantile":
                        change_range = (
                            protein_protein_interaction_network.
                            convert_quantile_to_change(
                                network,
                                subset.get("change",
                                           (0.025, 0.975))[0], subset["time"],
                                subset["post-translational modification"],
                                combination.SITE_COMBINATION[subset.get(
                                    "site combination", "absmax")]),
                            protein_protein_interaction_network.
                            convert_quantile_to_change(
                                network,
                                subset.get("change",
                                           (0.025, 0.975))[1], subset["time"],
                                subset["post-translational modification"],
                                combination.SITE_COMBINATION[subset.get(
                                    "site combination", "absmax")]))
                        proteins.intersection_update(
                            protein_protein_interaction_network.get_proteins(
                                network, subset["time"],
                                subset["post-translational modification"],
                                combination.SITE_COMBINATION[subset.get(
                                    "site combination", "absmax")],
                                lambda change: change <= change_range[
                                    0] or change >= change_range[1]))
                    else:
                        change_range = (subset.get("change", (-1.0, 1.0))[0],
                                        subset.get("change", (-1.0, 1.0))[1])
                        proteins.intersection_update(
                            protein_protein_interaction_network.get_proteins(
                                network, subset["time"],
                                subset["post-translational modification"],
                                combination.SITE_COMBINATION[subset.get(
                                    "site combination", "absmax")],
                                lambda change: change <= change_range[
                                    0] or change >= change_range[1]))

            reactome_network = pathway_network.get_pathway_network(
                nx.induced_subgraph(network, proteins),
                enrichment_test=test.ENRICHMENT_TEST[
                    configuration["Reactome network"].get(
                        "test", "hypergeometric")],
                multiple_testing_correction=correction.CORRECTION[
                    configuration["Reactome network"].get(
                        "correction", "Benjamini-Hochberg")],
                taxonomy_identifier=configuration["Reactome network"].get(
                    "taxonomy identifier", 9606))

        else:
            reactome_network = pathway_network.get_pathway_network(
                network,
                enrichment_test=test.ENRICHMENT_TEST[
                    configuration["Reactome network"].get(
                        "test", "hypergeometric")],
                multiple_testing_correction=correction.CORRECTION[
                    configuration["Reactome network"].get(
                        "correction", "Benjamini-Hochberg")],
                taxonomy_identifier=configuration["Reactome network"].get(
                    "taxonomy identifier", 9606))

        pathway_network.export(reactome_network,
                               "{}.reactome".format(logger.name),
                               ".{}".format(index) if index else "")

        if "Cytoscape" in configuration:
            styles.export(styles.get_pathway_network_style(reactome_network),
                          "{}.reactome".format(logger.name),
                          ".{}".format(index) if index else "")

    if "Gene Ontology network" in configuration:
        if "union" in configuration["Gene Ontology network"]:
            proteins = set()
            for subset in configuration["Gene Ontology network"]["union"]:
                if subset.get(
                        "time"
                ) in protein_protein_interaction_network.get_times(
                        network
                ) and subset.get(
                        "post-translational modification"
                ) in protein_protein_interaction_network.get_post_translational_modifications(
                        network, subset["time"]):
                    if subset.get("conversion") == "standard score":
                        change_range = (
                            protein_protein_interaction_network.
                            convert_standard_score_to_change(
                                network,
                                subset.get("change",
                                           (-2.0, 2.0))[0], subset["time"],
                                subset["post-translational modification"],
                                combination.SITE_COMBINATION[subset.get(
                                    "site combination", "absmax")]),
                            protein_protein_interaction_network.
                            convert_standard_score_to_change(
                                network,
                                subset.get("change",
                                           (-2.0, 2.0))[1], subset["time"],
                                subset["post-translational modification"],
                                combination.SITE_COMBINATION[subset.get(
                                    "site combination", "absmax")]))
                        proteins.update(
                            protein_protein_interaction_network.get_proteins(
                                network, subset["time"],
                                subset["post-translational modification"],
                                combination.SITE_COMBINATION[subset.get(
                                    "site combination", "absmax")],
                                lambda change: change <= change_range[
                                    0] or change >= change_range[1]))

                    elif subset.get("conversion") == "quantile":
                        change_range = (
                            protein_protein_interaction_network.
                            convert_quantile_to_change(
                                network,
                                subset.get("change",
                                           (0.025, 0.975))[0], subset["time"],
                                subset["post-translational modification"],
                                combination.SITE_COMBINATION[subset.get(
                                    "site combination", "absmax")]),
                            protein_protein_interaction_network.
                            convert_quantile_to_change(
                                network,
                                subset.get("change",
                                           (0.025, 0.975))[1], subset["time"],
                                subset["post-translational modification"],
                                combination.SITE_COMBINATION[subset.get(
                                    "site combination", "absmax")]))
                        proteins.update(
                            protein_protein_interaction_network.get_proteins(
                                network, subset["time"],
                                subset["post-translational modification"],
                                combination.SITE_COMBINATION[subset.get(
                                    "site combination", "absmax")],
                                lambda change: change <= change_range[
                                    0] or change >= change_range[1]))
                    else:
                        change_range = (subset.get("change", (-1.0, 1.0))[0],
                                        subset.get("change", (-1.0, 1.0))[1])
                        proteins.update(
                            protein_protein_interaction_network.get_proteins(
                                network, subset["time"],
                                subset["post-translational modification"],
                                combination.SITE_COMBINATION[subset.get(
                                    "site combination", "absmax")],
                                lambda change: change <= change_range[
                                    0] or change >= change_range[1]))

            gene_ontology_network = ontology_network.get_ontology_network(
                nx.induced_subgraph(network, proteins),
                namespaces=configuration["Gene Ontology network"].get(
                    "namespaces", [
                        "biological_process", "cellular_component",
                        "molecular_function"
                    ]),
                enrichment_test=test.ENRICHMENT_TEST[
                    configuration["Gene Ontology network"].get(
                        "test", "hypergeometric")],
                multiple_testing_correction=correction.CORRECTION[
                    configuration["Gene Ontology network"].get(
                        "correction", "Benjamini-Hochberg")],
                taxonomy_identifier=configuration["Gene Ontology network"].get(
                    "taxonomy identifier", 9606))

        elif "intersection" in configuration["Gene Ontology network"]:
            proteins = set(network.nodes())
            for subset in configuration["Gene Ontology network"][
                    "intersection"]:
                if subset.get(
                        "time"
                ) in protein_protein_interaction_network.get_times(
                        network
                ) and subset.get(
                        "post-translational modification"
                ) in protein_protein_interaction_network.get_post_translational_modifications(
                        network, subset["time"]):
                    if subset.get("conversion") == "standard score":
                        change_range = (
                            protein_protein_interaction_network.
                            convert_standard_score_to_change(
                                network,
                                subset.get("change",
                                           (-2.0, 2.0))[0], subset["time"],
                                subset["post-translational modification"],
                                combination.SITE_COMBINATION[subset.get(
                                    "site combination", "absmax")]),
                            protein_protein_interaction_network.
                            convert_standard_score_to_change(
                                network,
                                subset.get("change",
                                           (-2.0, 2.0))[1], subset["time"],
                                subset["post-translational modification"],
                                combination.SITE_COMBINATION[subset.get(
                                    "site combination", "absmax")]))
                        proteins.intersection_update(
                            protein_protein_interaction_network.get_proteins(
                                network, subset["time"],
                                subset["post-translational modification"],
                                combination.SITE_COMBINATION[subset.get(
                                    "site combination", "absmax")],
                                lambda change: change <= change_range[
                                    0] or change >= change_range[1]))

                    elif subset.get("conversion") == "quantile":
                        change_range = (
                            protein_protein_interaction_network.
                            convert_quantile_to_change(
                                network,
                                subset.get("change",
                                           (0.025, 0.975))[0], subset["time"],
                                subset["post-translational modification"],
                                combination.SITE_COMBINATION[subset.get(
                                    "site combination", "absmax")]),
                            protein_protein_interaction_network.
                            convert_quantile_to_change(
                                network,
                                subset.get("change",
                                           (0.025, 0.975))[1], subset["time"],
                                subset["post-translational modification"],
                                combination.SITE_COMBINATION[subset.get(
                                    "site combination", "absmax")]))
                        proteins.intersection_update(
                            protein_protein_interaction_network.get_proteins(
                                network, subset["time"],
                                subset["post-translational modification"],
                                combination.SITE_COMBINATION[subset.get(
                                    "site combination", "absmax")],
                                lambda change: change <= change_range[
                                    0] or change >= change_range[1]))
                    else:
                        change_range = (subset.get("change", (-1.0, 1.0))[0],
                                        subset.get("change", (-1.0, 1.0))[1])
                        proteins.intersection_update(
                            protein_protein_interaction_network.get_proteins(
                                network, subset["time"],
                                subset["post-translational modification"],
                                combination.SITE_COMBINATION[subset.get(
                                    "site combination", "absmax")],
                                lambda change: change <= change_range[
                                    0] or change >= change_range[1]))

            gene_ontology_network = ontology_network.get_ontology_network(
                nx.induced_subgraph(network, proteins),
                namespaces=configuration["Gene Ontology network"].get(
                    "namespaces", [
                        "biological_process", "cellular_component",
                        "molecular_function"
                    ]),
                enrichment_test=test.ENRICHMENT_TEST[
                    configuration["Gene Ontology network"].get(
                        "test", "hypergeometric")],
                multiple_testing_correction=correction.CORRECTION[
                    configuration["Gene Ontology network"].get(
                        "correction", "Benjamini-Hochberg")],
                taxonomy_identifier=configuration["Gene Ontology network"].get(
                    "taxonomy identifier", 9606))

        else:
            gene_ontology_network = ontology_network.get_ontology_network(
                network,
                namespaces=configuration["Gene Ontology network"].get(
                    "namespaces", [
                        "biological_process", "cellular_component",
                        "molecular_function"
                    ]),
                enrichment_test=test.ENRICHMENT_TEST[
                    configuration["Gene Ontology network"].get(
                        "test", "hypergeometric")],
                multiple_testing_correction=correction.CORRECTION[
                    configuration["Gene Ontology network"].get(
                        "correction", "Benjamini-Hochberg")],
                taxonomy_identifier=configuration["Gene Ontology network"].get(
                    "taxonomy identifier", 9606))

        ontology_network.export(gene_ontology_network,
                                "{}.gene_ontology".format(logger.name),
                                ".{}".format(index) if index else "")

        if "Cytoscape" in configuration:
            styles.export(
                styles.get_ontology_network_style(gene_ontology_network),
                "{}.go".format(logger.name),
                ".{}".format(index) if index else "")


def process_configuration_file(configuration_file: str) -> None:
    """
    Launches execution of a workflow.

    Args:
        configuration_file: file name of configuration file.
    """
    with open(configuration_file) as configuration:
        process_configuration(
            json.load(configuration),
            logging.getLogger(
                os.path.splitext(os.path.basename(configuration_file))[0]))


def main() -> None:
    parser = argparse.ArgumentParser()

    parser.add_argument("-c",
                        "--configurations",
                        help="configuration files",
                        nargs="+",
                        required=True)

    parser.add_argument(
        "-l",
        "--level",
        help="logging level (default: INFO)",
        type=str,
        default="INFO",
        choices=["CRITICAL", "ERROR", "WARNING", "INFO", "DEBUG"])

    parser.add_argument(
        "-p",
        "--processes",
        help="maximum number of concurrent processes (default: {})".format(
            os.cpu_count()),
        type=int,
        default=os.cpu_count())

    args = parser.parse_args()

    logging.basicConfig(stream=sys.stdout,
                        level=args.level,
                        format="%(name)s\t%(message)s")

    with concurrent.futures.ProcessPoolExecutor(
            max_workers=args.processes) as executor:
        for process in executor.map(process_configuration_file,
                                    args.configurations):
            process.results()


if __name__ == "__main__":
    main()
