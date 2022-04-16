"""pipeline"""
import argparse
import concurrent.futures
import json
import logging
import os
import re
import sys
from typing import Any, Optional

import networkx as nx

from cytoscape import (gene_ontology_network_style,
                       protein_protein_interaction_network_style,
                       reactome_network_style)
from interface import (combination, conversion, correction, default,
                       modularization, test)
from networks import (gene_ontology_network,
                      protein_protein_interaction_network, reactome_network)


def process_configuration(configurations: list[dict[str, Any]],
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
    network = protein_protein_interaction_network.get_network()

    for entry in configuration.get("genes", {}):
        if "file" in entry and "accession column" in entry:
            protein_protein_interaction_network.add_genes_from_table(
                network,
                file_name=entry["file"],
                gene_accession_column=entry["accession column"],
                gene_accession_format=re.compile(
                    entry.get("accession format", "^(.+?)$")),
                sheet_name=entry.get("sheet", 1) -
                1 if isinstance(entry.get("sheet", 1), int) else entry["sheet"],
                header=entry.get("header", 1) - 1,
                taxonomy_identifier=entry.get("species", 9606),
            )

        elif "accessions" in entry:
            network.add_nodes_from(entry["accessions"])

        protein_protein_interaction_network.map_genes(
            network, entry.get("species", 9606))

    for entry in configuration.get("proteins", {}):
        if "file" in entry and "accession column" in entry:
            protein_protein_interaction_network.add_proteins_from_table(
                network,
                file_name=entry["file"],
                protein_accession_column=entry["accession column"],
                protein_accession_format=re.compile(
                    entry.get("accession format", "^(.+?)$")),
                time=entry.get("time", 0),
                modification=entry.get("post-translational modification", "M"),
                position_column=entry.get("position column", ""),
                position_format=re.compile(
                    entry.get("position format", "^(.+?)$")),
                replicates=entry.get("replicate columns", []),
                sheet_name=entry.get("sheet", 1) -
                1 if isinstance(entry.get("sheet", 1), int) else entry["sheet"],
                header=entry.get("header", 1) - 1,
                number_sites=entry.get("sites", 5),
                number_replicates=entry.get("replicates", 1),
                replicate_combination=combination.REPLICATE_COMBINATION[
                    entry.get("replicate combination", "mean")],
                measurement_conversion=conversion.LOGARITHM[entry.get(
                    "logarithm")])

        elif "accessions" in entry:
            network.add_nodes_from(entry["accessions"])

    for entry in configuration.get("networks", []):
        network = nx.compose(network, nx.readwrite.graphml.read_graphml(entry))

    if configuration.get("proteins", {}) or configuration.get("networks", {}):
        protein_protein_interaction_network.map_proteins(network)

    if "protein-protein interactions" in configuration:
        for neighbors in range(
                max(configuration["protein-protein interactions"].get(
                    database, {}).get("neighbors", 0)
                    for database in ("BioGRID", "CORUM", "IntAct", "MINT",
                                     "Reactome", "STRING"))):
            proteins = set()
            if "BioGRID" in configuration[
                    "protein-protein interactions"] and configuration[
                        "protein-protein interactions"]["BioGRID"].get(
                            "neighbors", 0) > neighbors:
                proteins.update(
                    protein_protein_interaction_network.
                    get_neighbors_from_biogrid(
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
                                "species", 9606),
                        version=tuple(
                            int(v) for v in
                            configuration["protein-protein interactions"]
                            ["BioGRID"].get("version").split(".")) if "version"
                        in configuration["protein-protein interactions"]
                        ["BioGRID"] else None))

            if "CORUM" in configuration[
                    "protein-protein interactions"] and configuration[
                        "protein-protein interactions"]["CORUM"].get(
                            "neighbors", 0) > neighbors:
                proteins.update(
                    protein_protein_interaction_network.
                    get_neighbors_from_corum(
                        network,
                        purification_methods=configuration[
                            "protein-protein interactions"]["CORUM"].get(
                                "purification methods", [])))

            if "IntAct" in configuration[
                    "protein-protein interactions"] and configuration[
                        "protein-protein interactions"]["IntAct"].get(
                            "neighbors", 0) > neighbors:
                proteins.update(
                    protein_protein_interaction_network.
                    get_neighbors_from_intact(
                        network,
                        interaction_detection_methods=configuration[
                            "protein-protein interactions"]["IntAct"].get(
                                "interaction detection methods", []),
                        interaction_types=configuration[
                            "protein-protein interactions"]["IntAct"].get(
                                "interaction types", []),
                        psi_mi_score=configuration[
                            "protein-protein interactions"]["IntAct"].get(
                                "score", 0.0),
                        taxonomy_identifier=configuration[
                            "protein-protein interactions"]["IntAct"].get(
                                "species", 9606),
                    ))

            if "MINT" in configuration[
                    "protein-protein interactions"] and configuration[
                        "protein-protein interactions"]["MINT"].get(
                            "neighbors", 0) > neighbors:
                proteins.update(
                    protein_protein_interaction_network.get_neighbors_from_mint(
                        network,
                        interaction_detection_methods=configuration[
                            "protein-protein interactions"]["MINT"].get(
                                "interaction detection methods", []),
                        interaction_types=configuration[
                            "protein-protein interactions"]["MINT"].get(
                                "interaction types", []),
                        psi_mi_score=configuration[
                            "protein-protein interactions"]["MINT"].get(
                                "score", 0.0),
                        taxonomy_identifier=configuration[
                            "protein-protein interactions"]["MINT"].get(
                                "species", 9606),
                    ))

            if "Reactome" in configuration[
                    "protein-protein interactions"] and configuration[
                        "protein-protein interactions"]["Reactome"].get(
                            "neighbors", 0) > neighbors:
                proteins.update(
                    protein_protein_interaction_network.
                    get_neighbors_from_reactome(
                        network,
                        interaction_context=configuration[
                            "protein-protein interactions"]["Reactome"].get(
                                "interaction context", []),
                        interaction_type=configuration[
                            "protein-protein interactions"]["Reactome"].get(
                                "interaction type", []),
                        taxonomy_identifier=configuration[
                            "protein-protein interactions"]["Reactome"].get(
                                "species", 9606),
                    ))

            if "STRING" in configuration[
                    "protein-protein interactions"] and configuration[
                        "protein-protein interactions"]["STRING"].get(
                            "neighbors", 0) > neighbors:
                proteins.update(
                    protein_protein_interaction_network.
                    get_neighbors_from_string(
                        network,
                        neighborhood=configuration[
                            "protein-protein interactions"]["STRING"].get(
                                "neighborhood score", 0.0),
                        neighborhood_transferred=configuration[
                            "protein-protein interactions"]["STRING"].get(
                                "neighborhood transferred score", 0.0),
                        fusion=configuration["protein-protein interactions"]
                        ["STRING"].get("fusion score", 0.0),
                        cooccurence=configuration[
                            "protein-protein interactions"]["STRING"].get(
                                "cooccurence score", 0.0),
                        homology=configuration["protein-protein interactions"]
                        ["STRING"].get("homology score", 0.0),
                        coexpression=configuration[
                            "protein-protein interactions"]["STRING"].get(
                                "coexpression score", 0.0),
                        coexpression_transferred=configuration[
                            "protein-protein interactions"]["STRING"].get(
                                "coexpression transferred score", 0.0),
                        experiments=configuration[
                            "protein-protein interactions"]["STRING"].get(
                                "experiments score", 0.0),
                        experiments_transferred=configuration[
                            "protein-protein interactions"]["STRING"].get(
                                "experiments transferred score", 0.0),
                        database=configuration["protein-protein interactions"]
                        ["STRING"].get("database score", 0.0),
                        database_transferred=configuration[
                            "protein-protein interactions"]["STRING"].get(
                                "database transferred score", 0.0),
                        textmining=configuration["protein-protein interactions"]
                        ["STRING"].get("textmining score", 0.0),
                        textmining_transferred=configuration[
                            "protein-protein interactions"]["STRING"].get(
                                "textmining transferred score", 0.0),
                        combined_score=configuration[
                            "protein-protein interactions"]["STRING"].get(
                                "combined score", 0.0),
                        physical=configuration["protein-protein interactions"]
                        ["STRING"].get("physical", False),
                        taxonomy_identifier=configuration[
                            "protein-protein interactions"]["STRING"].get(
                                "species", 9606),
                        version=configuration["protein-protein interactions"]
                        ["STRING"].get("version", 11.5),
                    ))

            network.add_nodes_from(proteins)
            protein_protein_interaction_network.map_proteins(network)

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
                        "species", 9606),
                version=tuple(
                    int(v)
                    for v in configuration["protein-protein interactions"]
                    ["BioGRID"].get("version").split(".")) if "version"
                in configuration["protein-protein interactions"]["BioGRID"] else
                None)

        if "CORUM" in configuration["protein-protein interactions"]:
            protein_protein_interaction_network.add_protein_protein_interactions_from_corum(
                network,
                purification_methods=configuration[
                    "protein-protein interactions"]["CORUM"].get(
                        "purification methods", []))

        if "IntAct" in configuration["protein-protein interactions"]:
            protein_protein_interaction_network.add_protein_protein_interactions_from_intact(
                network,
                interaction_detection_methods=configuration[
                    "protein-protein interactions"]["IntAct"].get(
                        "interaction detection methods", []),
                interaction_types=configuration["protein-protein interactions"]
                ["IntAct"].get("interaction types", []),
                psi_mi_score=configuration["protein-protein interactions"]
                ["IntAct"].get("score", 0.0),
                taxonomy_identifier=configuration[
                    "protein-protein interactions"]["IntAct"].get(
                        "species", 9606),
            )

        if "MINT" in configuration["protein-protein interactions"]:
            protein_protein_interaction_network.add_protein_protein_interactions_from_mint(
                network,
                interaction_detection_methods=configuration[
                    "protein-protein interactions"]["MINT"].get(
                        "interaction detection methods", []),
                interaction_types=configuration["protein-protein interactions"]
                ["MINT"].get("interaction types", []),
                psi_mi_score=configuration["protein-protein interactions"]
                ["MINT"].get("score", 0.0),
                taxonomy_identifier=configuration[
                    "protein-protein interactions"]["MINT"].get(
                        "species", 9606),
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
                        "species", 9606),
            )

        if "STRING" in configuration["protein-protein interactions"]:
            protein_protein_interaction_network.add_protein_protein_interactions_from_string(
                network,
                neighborhood=configuration["protein-protein interactions"]
                ["STRING"].get("neighborhood score", 0.0),
                neighborhood_transferred=configuration[
                    "protein-protein interactions"]["STRING"].get(
                        "neighborhood transferred score", 0.0),
                fusion=configuration["protein-protein interactions"]
                ["STRING"].get("fusion score", 0.0),
                cooccurence=configuration["protein-protein interactions"]
                ["STRING"].get("cooccurence score", 0.0),
                homology=configuration["protein-protein interactions"]
                ["STRING"].get("homology score", 0.0),
                coexpression=configuration["protein-protein interactions"]
                ["STRING"].get("coexpression score", 0.0),
                coexpression_transferred=configuration[
                    "protein-protein interactions"]["STRING"].get(
                        "coexpression transferred score", 0.0),
                experiments=configuration["protein-protein interactions"]
                ["STRING"].get("experiments score", 0.0),
                experiments_transferred=configuration[
                    "protein-protein interactions"]["STRING"].get(
                        "experiments transferred score", 0.0),
                database=configuration["protein-protein interactions"]
                ["STRING"].get("database score", 0.0),
                database_transferred=configuration[
                    "protein-protein interactions"]["STRING"].get(
                        "database transferred score", 0.0),
                textmining=configuration["protein-protein interactions"]
                ["STRING"].get("textmining score", 0.0),
                textmining_transferred=configuration[
                    "protein-protein interactions"]["STRING"].get(
                        "textmining transferred score", 0.0),
                combined_score=configuration["protein-protein interactions"]
                ["STRING"].get("combined score", 0.0),
                physical=configuration["protein-protein interactions"]
                ["STRING"].get("physical", False),
                taxonomy_identifier=configuration[
                    "protein-protein interactions"]["STRING"].get(
                        "species", 9606),
                version=configuration["protein-protein interactions"]
                ["STRING"].get("version", 11.5),
            )

        if "Cytoscape" in configuration and any(
                protein_protein_interaction_network.
                get_post_translational_modifications(network, time)
                for time in protein_protein_interaction_network.get_times(
                    network)):
            style = protein_protein_interaction_network_style.get_style(
                network,
                bar_chart_range=default.MEASUREMENT_RANGE[
                    configuration["Cytoscape"].get("bar chart",
                                                   {}).get("conversion")],
                convert_measurement=conversion.MEASUREMENT_CONVERSION[
                    configuration["Cytoscape"].get("bar chart",
                                                   {}).get("conversion")],
                site_combination=combination.SITE_COMBINATION[
                    configuration["Cytoscape"].get("bar chart", {}).get(
                        "site combination", "maxabs")],
                confidence_score_combination=combination.
                CONFIDENCE_SCORE_COMBINATION[configuration["Cytoscape"].get(
                    "edge transparency")])

            protein_protein_interaction_network.set_edge_weights(
                network,
                weight=combination.CONFIDENCE_SCORE_COMBINATION[
                    configuration["Cytoscape"].get("edge transparency")],
                attribute="score")

            protein_protein_interaction_network_style.export(
                style,
                f"{logger.name}.{index}" if index else logger.name,
            )

            protein_protein_interaction_network.set_post_translational_modification(
                network)

            protein_protein_interaction_network.set_measurements(
                network,
                site_combination=combination.SITE_COMBINATION[
                    configuration["Cytoscape"].get("node color", {}).get(
                        "site combination", "maxabs")],
                measurements=default.MEASUREMENT_RANGE[
                    configuration["Cytoscape"].get("node color",
                                                   {}).get("conversion")],
                convert_measurement=conversion.MEASUREMENT_CONVERSION[
                    configuration["Cytoscape"].get("node color",
                                                   {}).get("conversion")],
            )

        protein_protein_interaction_network.export(
            network,
            f"{logger.name}.{index}" if index else logger.name,
        )

    if "Gene Ontology enrichment" in configuration:
        gene_ontology_enrichment = protein_protein_interaction_network.get_gene_ontology_enrichment(
            [network],
            enrichment_test=test.ENRICHMENT_TEST[
                configuration["Gene Ontology enrichment"].get(
                    "test", "hypergeometric")],
            multiple_testing_correction=correction.CORRECTION[
                configuration["Gene Ontology enrichment"].get(
                    "correction", "Benjamini-Hochberg")],
            taxonomy_identifier=configuration["Gene Ontology enrichment"].get(
                "species", 9606),
            namespaces=configuration["Gene Ontology enrichment"].get(
                "namespaces", [
                    "cellular_component", "molecular_function",
                    "biological_process"
                ]))

        for (term, name), p in sorted(gene_ontology_enrichment[network].items(),
                                      key=lambda item: item[1]):
            if p <= configuration["Gene Ontology enrichment"].get("p", 1.0):
                logger.info(f"{term}\t{p:.2e}\t{name}")

    if "Reactome enrichment" in configuration:
        reactome_enrichment = protein_protein_interaction_network.get_reactome_enrichment(
            [network],
            enrichment_test=test.ENRICHMENT_TEST[
                configuration["Reactome enrichment"].get(
                    "test", "hypergeometric")],
            multiple_testing_correction=correction.CORRECTION[
                configuration["Reactome enrichment"].get(
                    "correction", "Benjamini-Hochberg")],
            taxonomy_identifier=configuration["Reactome enrichment"].get(
                "species", 9606))

        for (pathway, name), p in sorted(reactome_enrichment[network].items(),
                                         key=lambda item: item[1]):
            if p <= configuration["Reactome enrichment"].get("p", 1.0):
                logger.info(f"{pathway}\t{p:.2e}\t{name}")

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
                    measurement_range = (
                        conversion.MEASUREMENT_CONVERSION[subset.get(
                            "conversion")]
                        (subset.get(
                            "measurement", default.MEASUREMENT_RANGE[subset.get(
                                "conversion")])[0],
                         protein_protein_interaction_network.get_measurements(
                             network, subset["time"],
                             subset["post-translational modification"],
                             combination.SITE_COMBINATION[subset.get(
                                 "site combination", "maxabs")])),
                        conversion.MEASUREMENT_CONVERSION[subset.get(
                            "conversion")]
                        (subset.get(
                            "measurement", default.MEASUREMENT_RANGE[subset.get(
                                "conversion")])[1],
                         protein_protein_interaction_network.get_measurements(
                             network, subset["time"],
                             subset["post-translational modification"],
                             combination.SITE_COMBINATION[subset.get(
                                 "site combination", "maxabs")])))

                    proteins.update(
                        protein_protein_interaction_network.get_proteins(
                            network, subset["time"],
                            subset["post-translational modification"],
                            combination.SITE_COMBINATION[subset.get(
                                "site combination", "maxabs")], lambda
                            measurement: measurement <= measurement_range[
                                0] or measurement >= measurement_range[1]))

            ontology_network = gene_ontology_network.get_network(
                proteins,
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
                    "species", 9606))

        elif "intersection" in configuration["Gene Ontology network"]:
            proteins = set(network)
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
                    measurement_range = (
                        conversion.MEASUREMENT_CONVERSION[subset.get(
                            "conversion")]
                        (subset.get(
                            "measurement", default.MEASUREMENT_RANGE[subset.get(
                                "conversion")])[0],
                         protein_protein_interaction_network.get_measurements(
                             network, subset["time"],
                             subset["post-translational modification"],
                             combination.SITE_COMBINATION[subset.get(
                                 "site combination", "maxabs")])),
                        conversion.MEASUREMENT_CONVERSION[subset.get(
                            "conversion")]
                        (subset.get(
                            "measurement", default.MEASUREMENT_RANGE[subset.get(
                                "conversion")])[1],
                         protein_protein_interaction_network.get_measurements(
                             network, subset["time"],
                             subset["post-translational modification"],
                             combination.SITE_COMBINATION[subset.get(
                                 "site combination", "maxabs")])))

                    proteins.intersection_update(
                        protein_protein_interaction_network.get_proteins(
                            network, subset["time"],
                            subset["post-translational modification"],
                            combination.SITE_COMBINATION[subset.get(
                                "site combination", "maxabs")], lambda
                            measurement: measurement <= measurement_range[
                                0] or measurement >= measurement_range[1]))

            ontology_network = gene_ontology_network.get_network(
                proteins,
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
                    "species", 9606))

        else:
            ontology_network = gene_ontology_network.get_network(
                network.nodes(),
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
                    "species", 9606))

        gene_ontology_network.export(
            ontology_network,
            f"{logger.name}.go.{index}" if index else f"{logger.name}.go")

        if "Cytoscape" in configuration:
            ontology_network_style = gene_ontology_network_style.get_style(
                ontology_network)
            gene_ontology_network_style.export(
                ontology_network_style,
                f"{logger.name}.go.{index}" if index else f"{logger.name}.go")

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
                    measurement_range = (
                        conversion.MEASUREMENT_CONVERSION[subset.get(
                            "conversion")]
                        (subset.get(
                            "measurement", default.MEASUREMENT_RANGE[subset.get(
                                "conversion")])[0],
                         protein_protein_interaction_network.get_measurements(
                             network, subset["time"],
                             subset["post-translational modification"],
                             combination.SITE_COMBINATION[subset.get(
                                 "site combination", "maxabs")])),
                        conversion.MEASUREMENT_CONVERSION[subset.get(
                            "conversion")]
                        (subset.get(
                            "measurement", default.MEASUREMENT_RANGE[subset.get(
                                "conversion")])[1],
                         protein_protein_interaction_network.get_measurements(
                             network, subset["time"],
                             subset["post-translational modification"],
                             combination.SITE_COMBINATION[subset.get(
                                 "site combination", "maxabs")])))

                    proteins.update(
                        protein_protein_interaction_network.get_proteins(
                            network, subset["time"],
                            subset["post-translational modification"],
                            combination.SITE_COMBINATION[subset.get(
                                "site combination", "maxabs")], lambda
                            measurement: measurement <= measurement_range[
                                0] or measurement >= measurement_range[1]))

            pathway_network = reactome_network.get_network(
                proteins,
                enrichment_test=test.ENRICHMENT_TEST[
                    configuration["Reactome network"].get(
                        "test", "hypergeometric")],
                multiple_testing_correction=correction.CORRECTION[
                    configuration["Reactome network"].get(
                        "correction", "Benjamini-Hochberg")],
                taxonomy_identifier=configuration["Reactome network"].get(
                    "species", 9606))

        elif "intersection" in configuration["Reactome network"]:
            proteins = set(network)
            for subset in configuration["Reactome network"]["intersection"]:
                if subset.get(
                        "time"
                ) in protein_protein_interaction_network.get_times(
                        network
                ) and subset.get(
                        "post-translational modification"
                ) in protein_protein_interaction_network.get_post_translational_modifications(
                        network, subset["time"]):
                    measurement_range = (
                        conversion.MEASUREMENT_CONVERSION[subset.get(
                            "conversion")]
                        (subset.get(
                            "measurement", default.MEASUREMENT_RANGE[subset.get(
                                "conversion")])[0],
                         protein_protein_interaction_network.get_measurements(
                             network, subset["time"],
                             subset["post-translational modification"],
                             combination.SITE_COMBINATION[subset.get(
                                 "site combination", "maxabs")])),
                        conversion.MEASUREMENT_CONVERSION[subset.get(
                            "conversion")]
                        (subset.get(
                            "measurement", default.MEASUREMENT_RANGE[subset.get(
                                "conversion")])[1],
                         protein_protein_interaction_network.get_measurements(
                             network, subset["time"],
                             subset["post-translational modification"],
                             combination.SITE_COMBINATION[subset.get(
                                 "site combination", "maxabs")])))

                    proteins.intersection_update(
                        protein_protein_interaction_network.get_proteins(
                            network, subset["time"],
                            subset["post-translational modification"],
                            combination.SITE_COMBINATION[subset.get(
                                "site combination", "maxabs")], lambda
                            measurement: measurement <= measurement_range[
                                0] or measurement >= measurement_range[1]))

            pathway_network = reactome_network.get_network(
                proteins,
                enrichment_test=test.ENRICHMENT_TEST[
                    configuration["Reactome network"].get(
                        "test", "hypergeometric")],
                multiple_testing_correction=correction.CORRECTION[
                    configuration["Reactome network"].get(
                        "correction", "Benjamini-Hochberg")],
                taxonomy_identifier=configuration["Reactome network"].get(
                    "species", 9606))

        else:
            pathway_network = reactome_network.get_network(
                network.nodes(),
                enrichment_test=test.ENRICHMENT_TEST[
                    configuration["Reactome network"].get(
                        "test", "hypergeometric")],
                multiple_testing_correction=correction.CORRECTION[
                    configuration["Reactome network"].get(
                        "correction", "Benjamini-Hochberg")],
                taxonomy_identifier=configuration["Reactome network"].get(
                    "species", 9606))

        reactome_network.export(
            pathway_network, f"{logger.name}.reactome.{index}"
            if index else f"{logger.name}.reactome")

        if "Cytoscape" in configuration:
            pathway_network_style = reactome_network_style.get_style(
                pathway_network)
            reactome_network_style.export(
                pathway_network_style, f"{logger.name}.reactome.{index}"
                if index else f"{logger.name}.reactome")

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
                            "species", 9606),
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
                            "species", 9606),
                        namespaces=configuration["module detection"]
                        ["Gene Ontology enrichment"]["network"].get(
                            "namespaces", [
                                "cellular_component", "molecular_function",
                                "biological_process"
                            ]),
                        annotation_as_reference=False)

        if "Reactome enrichment" in configuration["module detection"]:
            reactome_enrichment = {}
            if "annotation" in configuration["module detection"][
                    "Reactome enrichment"]:
                reactome_enrichment[
                    "annotation"] = protein_protein_interaction_network.get_reactome_enrichment(
                        modules,
                        enrichment_test=test.ENRICHMENT_TEST[
                            configuration["module detection"]
                            ["Reactome enrichment"]["annotation"].get(
                                "test", "hypergeometric")],
                        multiple_testing_correction=correction.CORRECTION[
                            configuration["module detection"]
                            ["Reactome enrichment"]["annotation"].get(
                                "correction", "Benjamini-Hochberg")],
                        taxonomy_identifier=configuration["module detection"]
                        ["Reactome enrichment"]["annotation"].get(
                            "species", 9606))

            if "network" in configuration["module detection"][
                    "Reactome enrichment"]:
                reactome_enrichment[
                    "network"] = protein_protein_interaction_network.get_reactome_enrichment(
                        modules,
                        enrichment_test=test.ENRICHMENT_TEST[
                            configuration["module detection"]
                            ["Reactome enrichment"]["network"].get(
                                "test", "hypergeometric")],
                        multiple_testing_correction=correction.CORRECTION[
                            configuration["module detection"]
                            ["Reactome enrichment"]["network"].get(
                                "correction", "Benjamini-Hochberg")],
                        taxonomy_identifier=configuration["module detection"]
                        ["Reactome enrichment"]["network"].get("species", 9606),
                        annotation_as_reference=False)

        if "measurement enrichment" in configuration["module detection"]:
            measurement_enrichment = {}
            if "proteins" in configuration["module detection"][
                    "measurement enrichment"]:
                measurement_enrichment[
                    "proteins"] = protein_protein_interaction_network.get_measurement_enrichment(
                        network,
                        modules,
                        measurements=default.MEASUREMENT_RANGE[
                            configuration["module detection"]
                            ["measurement enrichment"]["proteins"].get(
                                "conversion")],
                        convert_measurement=conversion.MEASUREMENT_CONVERSION[
                            configuration["module detection"]
                            ["measurement enrichment"]["proteins"].get(
                                "conversion")],
                        site_combination=combination.SITE_COMBINATION[
                            configuration["module detection"]
                            ["measurement enrichment"]["proteins"].get(
                                "site combination", "maxabs")],
                        enrichment_test=test.ENRICHMENT_TEST[
                            configuration["module detection"]
                            ["measurement enrichment"]["proteins"].get(
                                "test", "hypergeometric")],
                        multiple_testing_correction=correction.CORRECTION[
                            configuration["module detection"]
                            ["measurement enrichment"]["proteins"].get(
                                "correction", "Benjamini-Hochberg")])

            if "sites" in configuration["module detection"][
                    "measurement enrichment"]:
                measurement_enrichment[
                    "sites"] = protein_protein_interaction_network.get_measurement_enrichment(
                        network,
                        modules,
                        measurements=default.MEASUREMENT_RANGE[
                            configuration["module detection"]
                            ["measurement enrichment"]["sites"].get(
                                "conversion")],
                        convert_measurement=conversion.MEASUREMENT_CONVERSION[
                            configuration["module detection"]
                            ["measurement enrichment"]["sites"].get(
                                "conversion")],
                        enrichment_test=test.ENRICHMENT_TEST[
                            configuration["module detection"]
                            ["measurement enrichment"]["sites"].get(
                                "test", "hypergeometric")],
                        multiple_testing_correction=correction.CORRECTION[
                            configuration["module detection"]
                            ["measurement enrichment"]["sites"].get(
                                "correction", "Benjamini-Hochberg")])

        if "measurement location" in configuration["module detection"]:
            measurement_location = {}
            if "proteins" in configuration["module detection"][
                    "measurement location"]:
                measurement_location[
                    "proteins"] = protein_protein_interaction_network.get_measurement_location(
                        network,
                        modules,
                        combination.SITE_COMBINATION[
                            configuration["module detection"]
                            ["measurement location"]["proteins"].get(
                                "site combination", "maxabs")],
                        location_test=test.LOCATION_TEST[
                            configuration["module detection"]
                            ["measurement location"]["proteins"].get(
                                "test", "Wilcoxon")],
                        multiple_testing_correction=correction.CORRECTION[
                            configuration["module detection"]
                            ["measurement location"]["proteins"].get(
                                "correction", "Benjamini-Hochberg")])

            if "sites" in configuration["module detection"][
                    "measurement location"]:
                measurement_location[
                    "sites"] = protein_protein_interaction_network.get_measurement_location(
                        network,
                        modules,
                        location_test=test.LOCATION_TEST[
                            configuration["module detection"]
                            ["measurement location"]["sites"].get(
                                "test", "Wilcoxon")],
                        multiple_testing_correction=correction.CORRECTION[
                            configuration["module detection"]
                            ["measurement location"]["sites"].get(
                                "correction", "Benjamini-Hochberg")])

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
                            logger.info(
                                f"{k}\tannotation\t{term}\t{p:.2e}\t{name}")

                if "network" in configuration["module detection"][
                        "Gene Ontology enrichment"]:
                    for (term, name), p in sorted(
                            gene_ontology_enrichment["network"][module].items(),
                            key=lambda item: item[1]):
                        if p <= configuration["module detection"][
                                "Gene Ontology enrichment"]["network"].get(
                                    "p", 1.0):
                            export = True
                            logger.info(
                                f"{k}\tnetwork\t{term}\t{p:.2e}\t{name}")

            if "Reactome enrichment" in configuration["module detection"]:
                if "annotation" in configuration["module detection"][
                        "Reactome enrichment"]:
                    for (pathway, name), p in sorted(
                            reactome_enrichment["annotation"][module].items(),
                            key=lambda item: item[1]):
                        if p <= configuration["module detection"][
                                "Reactome enrichment"]["annotation"].get(
                                    "p", 1.0):
                            export = True
                            logger.info(
                                f"{k}\tannotation\t{pathway}\t{p:.2e}\t{name}")

                if "network" in configuration["module detection"][
                        "Reactome enrichment"]:
                    for (pathway, name), p in sorted(
                            reactome_enrichment["network"][module].items(),
                            key=lambda item: item[1]):
                        if p <= configuration["module detection"][
                                "Reactome enrichment"]["network"].get("p", 1.0):
                            export = True
                            logger.info(
                                f"{k}\tnetwork\t{pathway}\t{p:.2e}\t{name}")

            if "measurement enrichment" in configuration["module detection"]:
                if "proteins" in configuration["module detection"][
                        "measurement enrichment"]:
                    for time in measurement_enrichment["proteins"][module]:
                        for modification, p in sorted(
                                measurement_enrichment["proteins"][module]
                            [time].items(),
                                key=lambda item: item[1]):
                            if p <= configuration["module detection"][
                                    "measurement enrichment"]["proteins"].get(
                                        "p", 1.0):
                                export = True
                                logger.info(
                                    f"{k}\tprotein measurement enrichment\t"
                                    f"{time} {modification}\t{p:.2e}")

                if "sites" in configuration["module detection"][
                        "measurement enrichment"]:
                    for time in measurement_enrichment["sites"][module]:
                        for modification, p in sorted(
                                measurement_enrichment["sites"][module]
                            [time].items(),
                                key=lambda item: item[1]):
                            if p <= configuration["module detection"][
                                    "measurement enrichment"]["sites"].get(
                                        "p", 1.0):
                                export = True
                                logger.info(
                                    f"{k}\tsite measurement enrichment\t{time} "
                                    f"{modification}\t{p:.2e}")

            if "measurement location" in configuration["module detection"]:
                if "proteins" in configuration["module detection"][
                        "measurement location"]:
                    for time in measurement_location["proteins"][module]:
                        for modification, p in sorted(
                                measurement_location["proteins"][module]
                            [time].items(),
                                key=lambda item: item[1]):
                            if p <= configuration["module detection"][
                                    "measurement location"]["proteins"].get(
                                        "p", 1.0):
                                export = True
                                logger.info(
                                    f"{k}\tprotein measurement location\t"
                                    f"{time} {modification}\t{p:.2e}")

                if "sites" in configuration["module detection"][
                        "measurement location"]:
                    for time in measurement_location["sites"][module]:
                        for modification, p in sorted(
                                measurement_location["sites"][module]
                            [time].items(),
                                key=lambda item: item[1]):
                            if p <= configuration["module detection"][
                                    "measurement location"]["sites"].get(
                                        "p", 1.0):
                                export = True
                                logger.info(
                                    f"{k}\tsite measurement location\t{time} "
                                    f"{modification}\t{p:.2e}")

            if export:
                protein_protein_interaction_network.export(
                    module,
                    f"{logger.name}.{index}.{k}"
                    if index else f"{logger.name}.{k}",
                )


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
        help=
        f"maximum number of concurrent processes (default: {os.cpu_count()})",
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
            if process:
                process.results()


if __name__ == "__main__":
    main()
