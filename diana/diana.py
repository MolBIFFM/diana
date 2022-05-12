"""
DIANA: Data Integration And Network-based Analysis for post-translational
modification mass spectrometry data
"""
import argparse
import concurrent.futures
import json
import logging
import os
import re
import sys
from typing import Any, Mapping, Optional, Sequence

import networkx as nx

from cytoscape import (gene_ontology_network_style,
                       protein_interaction_network_style,
                       reactome_network_style)
from databases import corum, gene_ontology, reactome
from interface import (combination, conversion, correction, default,
                       modularization, test)
from networks import (gene_ontology_network, protein_interaction_network,
                      reactome_network)


def process_workflow(configuration: Mapping[str, Any],
                     logger: logging.Logger,
                     index: Optional[int | str] = None) -> None:
    """
    Executes a workflow specified in configuration.

    Args:
        configuration: The specification of a workflow.
        logger: A logger.
    """
    network = protein_interaction_network.get_network()

    for entry in configuration.get("genes", {}):
        if "file" in entry and "accession column" in entry:
            protein_interaction_network.add_genes_from_table(
                network,
                file_name=entry["file"],
                gene_accession_column=entry["accession column"],
                gene_accession_format=re.compile(
                    entry.get("accession format", "^(.+?)$")),
                sheet_name=entry.get("sheet", 1) -
                1 if isinstance(entry.get("sheet", 1), int) else entry["sheet"],
                header=entry.get("header", 1) - 1,
            )

        elif "accessions" in entry:
            network.add_nodes_from(entry["accessions"])

    for organism in set(
            entry.get("organism", 9606)
            for entry in configuration.get("genes", {})):
        protein_interaction_network.map_genes(network, organism)

    for entry in configuration.get("proteins", {}):
        if "file" in entry and "accession column" in entry:
            protein_interaction_network.add_proteins_from_table(
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

    for entry in configuration.get("networks", {}):
        network = nx.compose(
            network,
            nx.readwrite.graphml.read_graphml(entry["network"]) if
            entry.get("network") else protein_interaction_network.get_network())

    for organism in set(
            entry.get("organism", 9606)
            for entry in configuration.get("proteins", {})) | set(
                entry.get("organism", 9606)
                for entry in configuration.get("networks", {})):
        protein_interaction_network.map_proteins(network, organism)

    if "protein-protein interactions" in configuration:
        for neighbors in range(
                max(configuration["protein-protein interactions"].get(
                    database, {}).get("neighbors", 0) for database in
                    configuration["protein-protein interactions"])):
            interacting_proteins = set()
            if "BioGRID" in configuration[
                    "protein-protein interactions"] and configuration[
                        "protein-protein interactions"]["BioGRID"].get(
                            "neighbors", 0) > neighbors:
                interacting_proteins.update(
                    protein_interaction_network.get_neighbors_from_biogrid(
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
                        organism=configuration["protein-protein interactions"]
                        ["BioGRID"].get("organism", 9606),
                        version=configuration["protein-protein interactions"]
                        ["BioGRID"].get("version")))

            if "CORUM" in configuration[
                    "protein-protein interactions"] and configuration[
                        "protein-protein interactions"]["CORUM"].get(
                            "neighbors", 0) > neighbors:
                interacting_proteins.update(
                    protein_interaction_network.get_neighbors_from_corum(
                        network,
                        purification_methods=configuration[
                            "protein-protein interactions"]["CORUM"].get(
                                "purification methods", [])))

            if "IntAct" in configuration[
                    "protein-protein interactions"] and configuration[
                        "protein-protein interactions"]["IntAct"].get(
                            "neighbors", 0) > neighbors:
                interacting_proteins.update(
                    protein_interaction_network.get_neighbors_from_intact(
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
                        organism=configuration["protein-protein interactions"]
                        ["IntAct"].get("organism", 9606),
                    ))

            if "MINT" in configuration[
                    "protein-protein interactions"] and configuration[
                        "protein-protein interactions"]["MINT"].get(
                            "neighbors", 0) > neighbors:
                interacting_proteins.update(
                    protein_interaction_network.get_neighbors_from_mint(
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
                        organism=configuration["protein-protein interactions"]
                        ["MINT"].get("organism", 9606),
                    ))

            if "Reactome" in configuration[
                    "protein-protein interactions"] and configuration[
                        "protein-protein interactions"]["Reactome"].get(
                            "neighbors", 0) > neighbors:
                interacting_proteins.update(
                    protein_interaction_network.get_neighbors_from_reactome(
                        network,
                        interaction_context=configuration[
                            "protein-protein interactions"]["Reactome"].get(
                                "interaction context", []),
                        interaction_type=configuration[
                            "protein-protein interactions"]["Reactome"].get(
                                "interaction type", []),
                        organism=configuration["protein-protein interactions"]
                        ["Reactome"].get("organism", 9606),
                    ))

            if "STRING" in configuration[
                    "protein-protein interactions"] and configuration[
                        "protein-protein interactions"]["STRING"].get(
                            "neighbors", 0) > neighbors:
                interacting_proteins.update(
                    protein_interaction_network.get_neighbors_from_string(
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
                        organism=configuration["protein-protein interactions"]
                        ["STRING"].get("organism", 9606),
                        version=configuration["protein-protein interactions"]
                        ["STRING"].get("version", 11.5),
                    ))

            network.add_nodes_from(interacting_proteins)
            for organism in set(configuration["protein-protein interactions"]
                                [database].get("organism", 9606) for database in
                                configuration["protein-protein interactions"]):
                protein_interaction_network.map_proteins(network, organism)

        if "BioGRID" in configuration["protein-protein interactions"]:
            protein_interaction_network.add_protein_interactions_from_biogrid(
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
                organism=configuration["protein-protein interactions"]
                ["BioGRID"].get("organism", 9606),
                version=configuration["protein-protein interactions"]
                ["BioGRID"].get("version"))

        if "CORUM" in configuration["protein-protein interactions"]:
            protein_interaction_network.add_protein_interactions_from_corum(
                network,
                purification_methods=configuration[
                    "protein-protein interactions"]["CORUM"].get(
                        "purification methods", []))

        if "IntAct" in configuration["protein-protein interactions"]:
            protein_interaction_network.add_protein_interactions_from_intact(
                network,
                interaction_detection_methods=configuration[
                    "protein-protein interactions"]["IntAct"].get(
                        "interaction detection methods", []),
                interaction_types=configuration["protein-protein interactions"]
                ["IntAct"].get("interaction types", []),
                psi_mi_score=configuration["protein-protein interactions"]
                ["IntAct"].get("score", 0.0),
                organism=configuration["protein-protein interactions"]
                ["IntAct"].get("organism", 9606),
            )

        if "MINT" in configuration["protein-protein interactions"]:
            protein_interaction_network.add_protein_interactions_from_mint(
                network,
                interaction_detection_methods=configuration[
                    "protein-protein interactions"]["MINT"].get(
                        "interaction detection methods", []),
                interaction_types=configuration["protein-protein interactions"]
                ["MINT"].get("interaction types", []),
                psi_mi_score=configuration["protein-protein interactions"]
                ["MINT"].get("score", 0.0),
                organism=configuration["protein-protein interactions"]
                ["MINT"].get("organism", 9606),
            )

        if "Reactome" in configuration["protein-protein interactions"]:
            protein_interaction_network.add_protein_interactions_from_reactome(
                network,
                interaction_context=configuration[
                    "protein-protein interactions"]["Reactome"].get(
                        "interaction context", []),
                interaction_type=configuration["protein-protein interactions"]
                ["Reactome"].get("interaction type", []),
                organism=configuration["protein-protein interactions"]
                ["Reactome"].get("organism", 9606),
            )

        if "STRING" in configuration["protein-protein interactions"]:
            protein_interaction_network.add_protein_interactions_from_string(
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
                organism=configuration["protein-protein interactions"]
                ["STRING"].get("organism", 9606),
                version=configuration["protein-protein interactions"]
                ["STRING"].get("version", 11.5),
            )

        if "Cytoscape" in configuration and any(
                protein_interaction_network.get_modifications(network, time)
                for time in protein_interaction_network.get_times(network)):
            style = protein_interaction_network_style.get_style(
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

            protein_interaction_network.set_edge_weights(
                network,
                weight=combination.CONFIDENCE_SCORE_COMBINATION[
                    configuration["Cytoscape"].get("edge transparency")],
                attribute="score")

            protein_interaction_network_style.export(
                style,
                f"{logger.name}.{index}" if index else logger.name,
            )

            protein_interaction_network.set_post_translational_modification(
                network)

            protein_interaction_network.set_measurements(
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

        protein_interaction_network.export(
            network,
            f"{logger.name}.{index}" if index else logger.name,
        )

    if "CORUM enrichment" in configuration:
        if "subsets" in configuration["CORUM enrichment"]:
            if configuration["CORUM enrichment"].get("intersection", False):
                proteins = set(network.nodes())
            else:
                proteins = set()

            for subset in configuration["CORUM enrichment"]["subsets"]:
                if subset.get("time") in protein_interaction_network.get_times(
                        network) and subset.get(
                            "post-translational modification"
                        ) in protein_interaction_network.get_modifications(
                            network, subset["time"]):
                    measurement_range = (
                        conversion.MEASUREMENT_CONVERSION[subset.get(
                            "conversion")](
                                subset.get(
                                    "measurement", default.MEASUREMENT_RANGE[
                                        subset.get("conversion")])[0],
                                protein_interaction_network.get_measurements(
                                    network, subset["time"],
                                    subset["post-translational modification"],
                                    combination.SITE_COMBINATION[subset.get(
                                        "site combination", "maxabs")])),
                        conversion.MEASUREMENT_CONVERSION[subset.get(
                            "conversion")](
                                subset.get(
                                    "measurement", default.MEASUREMENT_RANGE[
                                        subset.get("conversion")])[1],
                                protein_interaction_network.get_measurements(
                                    network, subset["time"],
                                    subset["post-translational modification"],
                                    combination.SITE_COMBINATION[subset.get(
                                        "site combination", "maxabs")])))

                    if configuration["CORUM enrichment"].get(
                            "intersection", False):
                        proteins.intersection_update(
                            protein_interaction_network.get_proteins(
                                network,
                                subset["time"],
                                subset["post-translational modification"],
                                combination.SITE_COMBINATION[subset.get(
                                    "site combination", "maxabs")],
                                lambda measurement, measurement_range=
                                measurement_range: measurement <=
                                measurement_range[
                                    0] or measurement >= measurement_range[1]))

                    else:
                        proteins.update(
                            protein_interaction_network.get_proteins(
                                network,
                                subset["time"],
                                subset["post-translational modification"],
                                combination.SITE_COMBINATION[subset.get(
                                    "site combination", "maxabs")],
                                lambda measurement, measurement_range=
                                measurement_range: measurement <=
                                measurement_range[
                                    0] or measurement >= measurement_range[1]))

            corum_enrichment = corum.get_enrichment(
                [frozenset(proteins)],
                enrichment_test=test.ENRICHMENT_TEST[
                    configuration["CORUM enrichment"].get(
                        "test", "hypergeometric")],
                multiple_testing_correction=correction.CORRECTION[
                    configuration["CORUM enrichment"].get(
                        "correction", "Benjamini-Hochberg")],
                purification_methods=configuration["CORUM enrichment"].get(
                    "purification methods", []),
                organism=configuration["CORUM enrichment"].get(
                    "organism", 9606))

            for (protein_complex, name), p in sorted(
                    corum_enrichment[frozenset(proteins)].items(),
                    key=lambda item: item[1]):
                if p <= configuration["CORUM enrichment"].get("p", 1.0):
                    logger.info(f"{protein_complex}\t{p:.2e}\t{name}")

        else:
            corum_enrichment = corum.get_enrichment(
                [frozenset(network.nodes())],
                enrichment_test=test.ENRICHMENT_TEST[
                    configuration["CORUM enrichment"].get(
                        "test", "hypergeometric")],
                multiple_testing_correction=correction.CORRECTION[
                    configuration["CORUM enrichment"].get(
                        "correction", "Benjamini-Hochberg")],
                purification_methods=configuration["CORUM enrichment"].get(
                    "purification methods", []),
                organism=configuration["CORUM enrichment"].get(
                    "organism", 9606))

            for (protein_complex, name), p in sorted(corum_enrichment[frozenset(
                    network.nodes())].items(),
                                                     key=lambda item: item[1]):
                if p <= configuration["CORUM enrichment"].get("p", 1.0):
                    logger.info(f"{protein_complex}\t{p:.2e}\t{name}")

    if "Gene Ontology enrichment" in configuration:
        if "subsets" in configuration["Gene Ontology enrichment"]:
            if configuration["Gene Ontology enrichment"].get(
                    "intersection", False):
                proteins = set(network.nodes())
            else:
                set()

            for subset in configuration["Gene Ontology enrichment"]["subsets"]:
                if subset.get("time") in protein_interaction_network.get_times(
                        network) and subset.get(
                            "post-translational modification"
                        ) in protein_interaction_network.get_modifications(
                            network, subset["time"]):
                    measurement_range = (
                        conversion.MEASUREMENT_CONVERSION[subset.get(
                            "conversion")](
                                subset.get(
                                    "measurement", default.MEASUREMENT_RANGE[
                                        subset.get("conversion")])[0],
                                protein_interaction_network.get_measurements(
                                    network, subset["time"],
                                    subset["post-translational modification"],
                                    combination.SITE_COMBINATION[subset.get(
                                        "site combination", "maxabs")])),
                        conversion.MEASUREMENT_CONVERSION[subset.get(
                            "conversion")](
                                subset.get(
                                    "measurement", default.MEASUREMENT_RANGE[
                                        subset.get("conversion")])[1],
                                protein_interaction_network.get_measurements(
                                    network, subset["time"],
                                    subset["post-translational modification"],
                                    combination.SITE_COMBINATION[subset.get(
                                        "site combination", "maxabs")])))

                    if configuration["Gene Ontology enrichment"].get(
                            "intersection", False):
                        proteins.intersection_update(
                            protein_interaction_network.get_proteins(
                                network,
                                subset["time"],
                                subset["post-translational modification"],
                                combination.SITE_COMBINATION[subset.get(
                                    "site combination", "maxabs")],
                                lambda measurement, measurement_range=
                                measurement_range: measurement <=
                                measurement_range[
                                    0] or measurement >= measurement_range[1]))

                    else:
                        proteins.update(
                            protein_interaction_network.get_proteins(
                                network,
                                subset["time"],
                                subset["post-translational modification"],
                                combination.SITE_COMBINATION[subset.get(
                                    "site combination", "maxabs")],
                                lambda measurement, measurement_range=
                                measurement_range: measurement <=
                                measurement_range[
                                    0] or measurement >= measurement_range[1]))

            gene_ontology_enrichment = gene_ontology.get_enrichment(
                [frozenset(network.nodes())],
                enrichment_test=test.ENRICHMENT_TEST[
                    configuration["Gene Ontology enrichment"].get(
                        "test", "hypergeometric")],
                multiple_testing_correction=correction.CORRECTION[
                    configuration["Gene Ontology enrichment"].get(
                        "correction", "Benjamini-Hochberg")],
                organism=configuration["Gene Ontology enrichment"].get(
                    "organism", 9606),
                namespaces=configuration["Gene Ontology enrichment"].get(
                    "namespaces", [
                        "cellular_component", "molecular_function",
                        "biological_process"
                    ]))

            for (term, name), p in sorted(
                    gene_ontology_enrichment[frozenset(proteins)].items(),
                    key=lambda item: item[1]):
                if p <= configuration["Gene Ontology enrichment"].get("p", 1.0):
                    logger.info(f"{term}\t{p:.2e}\t{name}")

        else:
            gene_ontology_enrichment = gene_ontology.get_enrichment(
                [frozenset(network.nodes())],
                enrichment_test=test.ENRICHMENT_TEST[
                    configuration["Gene Ontology enrichment"].get(
                        "test", "hypergeometric")],
                multiple_testing_correction=correction.CORRECTION[
                    configuration["Gene Ontology enrichment"].get(
                        "correction", "Benjamini-Hochberg")],
                organism=configuration["Gene Ontology enrichment"].get(
                    "organism", 9606),
                namespaces=configuration["Gene Ontology enrichment"].get(
                    "namespaces", [
                        "cellular_component", "molecular_function",
                        "biological_process"
                    ]))

            for (term, name), p in sorted(gene_ontology_enrichment[frozenset(
                    network.nodes())].items(),
                                          key=lambda item: item[1]):
                if p <= configuration["Gene Ontology enrichment"].get("p", 1.0):
                    logger.info(f"{term}\t{p:.2e}\t{name}")

    if "Reactome enrichment" in configuration:
        if "subsets" in configuration["Reactome enrichment"]:
            if configuration["Reactome enrichment"].get("intersection", False):
                proteins = set(network.nodes())
            else:
                proteins = set()

            for subset in configuration["Reactome enrichment"]["subsets"]:
                if subset.get("time") in protein_interaction_network.get_times(
                        network) and subset.get(
                            "post-translational modification"
                        ) in protein_interaction_network.get_modifications(
                            network, subset["time"]):
                    measurement_range = (
                        conversion.MEASUREMENT_CONVERSION[subset.get(
                            "conversion")](
                                subset.get(
                                    "measurement", default.MEASUREMENT_RANGE[
                                        subset.get("conversion")])[0],
                                protein_interaction_network.get_measurements(
                                    network, subset["time"],
                                    subset["post-translational modification"],
                                    combination.SITE_COMBINATION[subset.get(
                                        "site combination", "maxabs")])),
                        conversion.MEASUREMENT_CONVERSION[subset.get(
                            "conversion")](
                                subset.get(
                                    "measurement", default.MEASUREMENT_RANGE[
                                        subset.get("conversion")])[1],
                                protein_interaction_network.get_measurements(
                                    network, subset["time"],
                                    subset["post-translational modification"],
                                    combination.SITE_COMBINATION[subset.get(
                                        "site combination", "maxabs")])))

                    if configuration["Reactome enrichment"].get(
                            "intersection", False):
                        proteins.intersection_update(
                            protein_interaction_network.get_proteins(
                                network,
                                subset["time"],
                                subset["post-translational modification"],
                                combination.SITE_COMBINATION[subset.get(
                                    "site combination", "maxabs")],
                                lambda measurement, measurement_range=
                                measurement_range: measurement <=
                                measurement_range[
                                    0] or measurement >= measurement_range[1]))

                    else:
                        proteins.update(
                            protein_interaction_network.get_proteins(
                                network,
                                subset["time"],
                                subset["post-translational modification"],
                                combination.SITE_COMBINATION[subset.get(
                                    "site combination", "maxabs")],
                                lambda measurement, measurement_range=
                                measurement_range: measurement <=
                                measurement_range[
                                    0] or measurement >= measurement_range[1]))

            reactome_enrichment = reactome.get_enrichment(
                [frozenset(network.nodes())],
                enrichment_test=test.ENRICHMENT_TEST[
                    configuration["Reactome enrichment"].get(
                        "test", "hypergeometric")],
                multiple_testing_correction=correction.CORRECTION[
                    configuration["Reactome enrichment"].get(
                        "correction", "Benjamini-Hochberg")],
                organism=configuration["Reactome enrichment"].get(
                    "organism", 9606))

            for (pathway, name), p in sorted(
                    reactome_enrichment[frozenset(proteins)].items(),
                    key=lambda item: item[1]):
                if p <= configuration["Reactome enrichment"].get("p", 1.0):
                    logger.info(f"{pathway}\t{p:.2e}\t{name}")
        else:
            reactome_enrichment = reactome.get_enrichment(
                [frozenset(network.nodes())],
                enrichment_test=test.ENRICHMENT_TEST[
                    configuration["Reactome enrichment"].get(
                        "test", "hypergeometric")],
                multiple_testing_correction=correction.CORRECTION[
                    configuration["Reactome enrichment"].get(
                        "correction", "Benjamini-Hochberg")],
                organism=configuration["Reactome enrichment"].get(
                    "organism", 9606))

            for (pathway, name), p in sorted(reactome_enrichment[frozenset(
                    network.nodes())].items(),
                                             key=lambda item: item[1]):
                if p <= configuration["Reactome enrichment"].get("p", 1.0):
                    logger.info(f"{pathway}\t{p:.2e}\t{name}")

    if "Gene Ontology network" in configuration:
        if "subsets" in configuration["Gene Ontology network"]:
            proteins = set(
                network.nodes()) if configuration["Gene Ontology network"].get(
                    "intersection", False) else set()
            for subset in configuration["Gene Ontology network"]["subsets"]:
                if subset.get("time") in protein_interaction_network.get_times(
                        network) and subset.get(
                            "post-translational modification"
                        ) in protein_interaction_network.get_modifications(
                            network, subset["time"]):
                    measurement_range = (
                        conversion.MEASUREMENT_CONVERSION[subset.get(
                            "conversion")](
                                subset.get(
                                    "measurement", default.MEASUREMENT_RANGE[
                                        subset.get("conversion")])[0],
                                protein_interaction_network.get_measurements(
                                    network, subset["time"],
                                    subset["post-translational modification"],
                                    combination.SITE_COMBINATION[subset.get(
                                        "site combination", "maxabs")])),
                        conversion.MEASUREMENT_CONVERSION[subset.get(
                            "conversion")](
                                subset.get(
                                    "measurement", default.MEASUREMENT_RANGE[
                                        subset.get("conversion")])[1],
                                protein_interaction_network.get_measurements(
                                    network, subset["time"],
                                    subset["post-translational modification"],
                                    combination.SITE_COMBINATION[subset.get(
                                        "site combination", "maxabs")])))

                    if configuration["Gene Ontology network"].get(
                            "intersection", False):
                        proteins.intersection_update(
                            protein_interaction_network.get_proteins(
                                network,
                                subset["time"],
                                subset["post-translational modification"],
                                combination.SITE_COMBINATION[subset.get(
                                    "site combination", "maxabs")],
                                lambda measurement, measurement_range=
                                measurement_range: measurement <=
                                measurement_range[
                                    0] or measurement >= measurement_range[1]))

                    else:
                        proteins.update(
                            protein_interaction_network.get_proteins(
                                network,
                                subset["time"],
                                subset["post-translational modification"],
                                combination.SITE_COMBINATION[subset.get(
                                    "site combination", "maxabs")],
                                lambda measurement, measurement_range=
                                measurement_range: measurement <=
                                measurement_range[
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
                organism=configuration["Gene Ontology network"].get(
                    "organism", 9606))

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
                organism=configuration["Gene Ontology network"].get(
                    "organism", 9606))

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
        if "subsets" in configuration["Reactome network"]:
            proteins = set(
                network.nodes()) if configuration["Reactome network"].get(
                    "intersection", False) else set()
            for subset in configuration["Reactome network"]["subsets"]:
                if subset.get("time") in protein_interaction_network.get_times(
                        network) and subset.get(
                            "post-translational modification"
                        ) in protein_interaction_network.get_modifications(
                            network, subset["time"]):
                    measurement_range = (
                        conversion.MEASUREMENT_CONVERSION[subset.get(
                            "conversion")](
                                subset.get(
                                    "measurement", default.MEASUREMENT_RANGE[
                                        subset.get("conversion")])[0],
                                protein_interaction_network.get_measurements(
                                    network, subset["time"],
                                    subset["post-translational modification"],
                                    combination.SITE_COMBINATION[subset.get(
                                        "site combination", "maxabs")])),
                        conversion.MEASUREMENT_CONVERSION[subset.get(
                            "conversion")](
                                subset.get(
                                    "measurement", default.MEASUREMENT_RANGE[
                                        subset.get("conversion")])[1],
                                protein_interaction_network.get_measurements(
                                    network, subset["time"],
                                    subset["post-translational modification"],
                                    combination.SITE_COMBINATION[subset.get(
                                        "site combination", "maxabs")])))

                    if configuration["Reactome network"].get(
                            "intersection", False):
                        proteins.intersection_update(
                            protein_interaction_network.get_proteins(
                                network,
                                subset["time"],
                                subset["post-translational modification"],
                                combination.SITE_COMBINATION[subset.get(
                                    "site combination", "maxabs")],
                                lambda measurement, measurement_range=
                                measurement_range: measurement <=
                                measurement_range[
                                    0] or measurement >= measurement_range[1]))

                    else:
                        proteins.update(
                            protein_interaction_network.get_proteins(
                                network,
                                subset["time"],
                                subset["post-translational modification"],
                                combination.SITE_COMBINATION[subset.get(
                                    "site combination", "maxabs")],
                                lambda measurement, measurement_range=
                                measurement_range: measurement <=
                                measurement_range[
                                    0] or measurement >= measurement_range[1]))

            pathway_network = reactome_network.get_network(
                proteins,
                enrichment_test=test.ENRICHMENT_TEST[
                    configuration["Reactome network"].get(
                        "test", "hypergeometric")],
                multiple_testing_correction=correction.CORRECTION[
                    configuration["Reactome network"].get(
                        "correction", "Benjamini-Hochberg")],
                organism=configuration["Reactome network"].get(
                    "organism", 9606))

        else:
            pathway_network = reactome_network.get_network(
                network.nodes(),
                enrichment_test=test.ENRICHMENT_TEST[
                    configuration["Reactome network"].get(
                        "test", "hypergeometric")],
                multiple_testing_correction=correction.CORRECTION[
                    configuration["Reactome network"].get(
                        "correction", "Benjamini-Hochberg")],
                organism=configuration["Reactome network"].get(
                    "organism", 9606))

        reactome_network.export(
            pathway_network, f"{logger.name}.reactome.{index}"
            if index else f"{logger.name}.reactome")

        if "Cytoscape" in configuration:
            pathway_network_style = reactome_network_style.get_style(
                pathway_network)
            reactome_network_style.export(
                pathway_network_style, f"{logger.name}.reactome.{index}"
                if index else f"{logger.name}.reactome")

    if "community detection" in configuration:
        protein_interaction_network.set_edge_weights(
            network,
            weight=combination.CONFIDENCE_SCORE_COMBINATION[
                configuration["community detection"].get("edge weight")])

        communities = protein_interaction_network.get_communities(
            network,
            community_size=configuration["community detection"].get(
                "community size", network.number_of_nodes()),
            community_size_combination=combination.MODULE_SIZE_COMBINATION[
                configuration["community detection"].get(
                    "community size combination", "mean")],
            algorithm=modularization.ALGORITHM[
                configuration["community detection"].get(
                    "algorithm", "Louvain")],
            resolution=configuration["community detection"].get(
                "resolution", 1.0))

        protein_interaction_network.remove_edge_weights(network)

        export = {community: False for community in communities}

        if "CORUM enrichment" in configuration["community detection"]:
            if "subsets" in configuration["community detection"][
                    "CORUM enrichment"]:
                if configuration["community detection"]["CORUM enrichment"].get(
                        "intersection", False):
                    cmty_proteins = {
                        community: set(community.nodes())
                        for community in communities
                    }
                else:
                    cmty_proteins = {
                        community: set() for community in communities
                    }

                for subset in configuration["community detection"][
                        "CORUM enrichment"]["subsets"]:
                    if subset.get(
                            "time") in protein_interaction_network.get_times(
                                network
                            ) and subset.get(
                                "post-translational modification"
                            ) in protein_interaction_network.get_modifications(
                                network, subset["time"]):
                        for community in cmty_proteins:
                            measurement_range = (
                                conversion.MEASUREMENT_CONVERSION[subset.get(
                                    "conversion")]
                                (subset.get(
                                    "measurement", default.MEASUREMENT_RANGE[
                                        subset.get("conversion")])[0],
                                 protein_interaction_network.get_measurements(
                                     community, subset["time"],
                                     subset["post-translational modification"],
                                     combination.SITE_COMBINATION[subset.get(
                                         "site combination", "maxabs")])),
                                conversion.MEASUREMENT_CONVERSION[subset.get(
                                    "conversion")]
                                (subset.get(
                                    "measurement", default.MEASUREMENT_RANGE[
                                        subset.get("conversion")])[1],
                                 protein_interaction_network.get_measurements(
                                     community, subset["time"],
                                     subset["post-translational modification"],
                                     combination.SITE_COMBINATION[subset.get(
                                         "site combination", "maxabs")])))

                            if configuration["community detection"][
                                    "CORUM enrichment"].get(
                                        "intersection", False):
                                cmty_proteins[community].intersection_update(
                                    protein_interaction_network.get_proteins(
                                        community,
                                        subset["time"],
                                        subset[
                                            "post-translational modification"],
                                        combination.SITE_COMBINATION[subset.get(
                                            "site combination", "maxabs")],
                                        lambda measurement, measurement_range=
                                        measurement_range: measurement <=
                                        measurement_range[0] or measurement >=
                                        measurement_range[1]))

                            else:
                                cmty_proteins[community].update(
                                    protein_interaction_network.get_proteins(
                                        community,
                                        subset["time"],
                                        subset[
                                            "post-translational modification"],
                                        combination.SITE_COMBINATION[subset.get(
                                            "site combination", "maxabs")],
                                        lambda measurement, measurement_range=
                                        measurement_range: measurement <=
                                        measurement_range[0] or measurement >=
                                        measurement_range[1]))

                corum_enrichment = corum.get_enrichment(
                    [
                        frozenset(cmty_proteins[community])
                        for community in communities
                    ],
                    enrichment_test=test.ENRICHMENT_TEST[
                        configuration["CORUM enrichment"].get(
                            "test", "hypergeometric")],
                    multiple_testing_correction=correction.CORRECTION[
                        configuration["CORUM enrichment"].get(
                            "correction", "Benjamini-Hochberg")],
                    purification_methods=configuration["CORUM enrichment"].get(
                        "purification methods", []),
                    organism=configuration["CORUM enrichment"].get(
                        "organism", 9606),
                    annotation_as_reference=configuration["community detection"]
                    ["CORUM enrichment"].get("annotation", False))

                for k, community in enumerate(sorted(
                        communities,
                        key=lambda community: community.number_of_nodes(),
                        reverse=True),
                                              start=1):
                    for (protein_complex,
                         name), p in sorted(corum_enrichment[frozenset(
                             cmty_proteins[community])].items(),
                                            key=lambda item: item[1]):
                        if p <= configuration["community detection"][
                                "CORUM enrichment"].get("p", 1.0):
                            export[community] = True
                            logger.info(f"{k}\t{protein_complex}\t{p:.2e}\t"
                                        f"{name}")
            else:
                corum_enrichment = corum.get_enrichment(
                    [frozenset(community.nodes()) for community in communities],
                    enrichment_test=test.ENRICHMENT_TEST[
                        configuration["community detection"]
                        ["CORUM enrichment"].get("test", "hypergeometric")],
                    multiple_testing_correction=correction.CORRECTION[
                        configuration["community detection"]
                        ["CORUM enrichment"].get("correction",
                                                 "Benjamini-Hochberg")],
                    purification_methods=configuration["community detection"]
                    ["CORUM enrichment"].get("purification methods", []),
                    organism=configuration["community detection"]
                    ["CORUM enrichment"].get("organism", 9606),
                    annotation_as_reference=configuration["community detection"]
                    ["CORUM enrichment"].get("annotation", False))

                for k, community in enumerate(sorted(
                        communities,
                        key=lambda community: community.number_of_nodes(),
                        reverse=True),
                                              start=1):
                    for (protein_complex,
                         name), p in sorted(corum_enrichment[frozenset(
                             community.nodes())].items(),
                                            key=lambda item: item[1]):
                        if p <= configuration["community detection"][
                                "CORUM enrichment"].get("p", 1.0):
                            export[community] = True
                            logger.info(f"{k}\t{protein_complex}\t{p:.2e}\t"
                                        f"{name}")

        if "Gene Ontology enrichment" in configuration["community detection"]:
            if "subsets" in configuration["community detection"][
                    "Gene Ontology enrichment"]:
                if configuration["community detection"][
                        "Gene Ontology enrichment"].get("intersection", False):
                    cmty_proteins = {
                        community: set(community.nodes())
                        for community in communities
                    }
                else:
                    cmty_proteins = {
                        community: set() for community in communities
                    }

                for subset in configuration["community detection"][
                        "Gene Ontology enrichment"]["subsets"]:
                    if subset.get(
                            "time") in protein_interaction_network.get_times(
                                network
                            ) and subset.get(
                                "post-translational modification"
                            ) in protein_interaction_network.get_modifications(
                                network, subset["time"]):
                        for community in cmty_proteins:
                            measurement_range = (
                                conversion.MEASUREMENT_CONVERSION[subset.get(
                                    "conversion")]
                                (subset.get(
                                    "measurement", default.MEASUREMENT_RANGE[
                                        subset.get("conversion")])[0],
                                 protein_interaction_network.get_measurements(
                                     community, subset["time"],
                                     subset["post-translational modification"],
                                     combination.SITE_COMBINATION[subset.get(
                                         "site combination", "maxabs")])),
                                conversion.MEASUREMENT_CONVERSION[subset.get(
                                    "conversion")]
                                (subset.get(
                                    "measurement", default.MEASUREMENT_RANGE[
                                        subset.get("conversion")])[1],
                                 protein_interaction_network.get_measurements(
                                     community, subset["time"],
                                     subset["post-translational modification"],
                                     combination.SITE_COMBINATION[subset.get(
                                         "site combination", "maxabs")])))

                            if configuration["community detection"][
                                    "Gene Ontology enrichment"].get(
                                        "intersection", False):
                                cmty_proteins[community].intersection_update(
                                    protein_interaction_network.get_proteins(
                                        community,
                                        subset["time"],
                                        subset[
                                            "post-translational modification"],
                                        combination.SITE_COMBINATION[subset.get(
                                            "site combination", "maxabs")],
                                        lambda measurement, measurement_range=
                                        measurement_range: measurement <=
                                        measurement_range[0] or measurement >=
                                        measurement_range[1]))

                            else:
                                cmty_proteins[community].update(
                                    protein_interaction_network.get_proteins(
                                        community,
                                        subset["time"],
                                        subset[
                                            "post-translational modification"],
                                        combination.SITE_COMBINATION[subset.get(
                                            "site combination", "maxabs")],
                                        lambda measurement, measurement_range=
                                        measurement_range: measurement <=
                                        measurement_range[0] or measurement >=
                                        measurement_range[1]))

                gene_ontology_enrichment = gene_ontology.get_enrichment(
                    [
                        frozenset(cmty_proteins[community])
                        for community in communities
                    ],
                    enrichment_test=test.ENRICHMENT_TEST[configuration[
                        "community detection"]["Gene Ontology enrichment"].get(
                            "test", "hypergeometric")],
                    multiple_testing_correction=correction.CORRECTION[
                        configuration["community detection"]
                        ["Gene Ontology enrichment"].get(
                            "correction", "Benjamini-Hochberg")],
                    organism=configuration["community detection"]
                    ["Gene Ontology enrichment"].get("organism", 9606),
                    namespaces=configuration["community detection"]
                    ["Gene Ontology enrichment"].get("namespaces", [
                        "cellular_component", "molecular_function",
                        "biological_process"
                    ]),
                    annotation_as_reference=configuration["community detection"]
                    ["Gene Ontology enrichment"].get("annotation", False))

                for k, community in enumerate(sorted(
                        communities,
                        key=lambda community: community.number_of_nodes(),
                        reverse=True),
                                              start=1):
                    for (term,
                         name), p in sorted(gene_ontology_enrichment[frozenset(
                             cmty_proteins[community])].items(),
                                            key=lambda item: item[1]):
                        if p <= configuration["community detection"][
                                "Gene Ontology enrichment"].get("p", 1.0):
                            export[community] = True
                            logger.info(f"{k}\t{term}\t{p:.2e}\t{name}")

            else:
                gene_ontology_enrichment = gene_ontology.get_enrichment(
                    [frozenset(community.nodes()) for community in communities],
                    enrichment_test=test.ENRICHMENT_TEST[configuration[
                        "community detection"]["Gene Ontology enrichment"].get(
                            "test", "hypergeometric")],
                    multiple_testing_correction=correction.CORRECTION[
                        configuration["community detection"]
                        ["Gene Ontology enrichment"].get(
                            "correction", "Benjamini-Hochberg")],
                    organism=configuration["community detection"]
                    ["Gene Ontology enrichment"].get("organism", 9606),
                    namespaces=configuration["community detection"]
                    ["Gene Ontology enrichment"].get("namespaces", [
                        "cellular_component", "molecular_function",
                        "biological_process"
                    ]),
                    annotation_as_reference=configuration["community detection"]
                    ["Gene Ontology enrichment"].get("annotation", False))

                for k, community in enumerate(sorted(
                        communities,
                        key=lambda community: community.number_of_nodes(),
                        reverse=True),
                                              start=1):
                    for (term,
                         name), p in sorted(gene_ontology_enrichment[frozenset(
                             community.nodes())].items(),
                                            key=lambda item: item[1]):
                        if p <= configuration["community detection"][
                                "Gene Ontology enrichment"].get("p", 1.0):
                            export[community] = True
                            logger.info(f"{k}\t{term}\t{p:.2e}\t{name}")

        if "Reactome enrichment" in configuration["community detection"]:
            if "subsets" in configuration["community detection"][
                    "Reactome enrichment"]:
                if configuration["community detection"][
                        "Reactome enrichment"].get("intersection", False):
                    cmty_proteins = {
                        community: set(community.nodes())
                        for community in communities
                    }
                else:
                    cmty_proteins = {
                        community: set() for community in communities
                    }

                for subset in configuration["community detection"][
                        "Reactome enrichment"]["subsets"]:
                    if subset.get(
                            "time") in protein_interaction_network.get_times(
                                network
                            ) and subset.get(
                                "post-translational modification"
                            ) in protein_interaction_network.get_modifications(
                                network, subset["time"]):
                        for community in cmty_proteins:
                            measurement_range = (
                                conversion.MEASUREMENT_CONVERSION[subset.get(
                                    "conversion")]
                                (subset.get(
                                    "measurement", default.MEASUREMENT_RANGE[
                                        subset.get("conversion")])[0],
                                 protein_interaction_network.get_measurements(
                                     community, subset["time"],
                                     subset["post-translational modification"],
                                     combination.SITE_COMBINATION[subset.get(
                                         "site combination", "maxabs")])),
                                conversion.MEASUREMENT_CONVERSION[subset.get(
                                    "conversion")]
                                (subset.get(
                                    "measurement", default.MEASUREMENT_RANGE[
                                        subset.get("conversion")])[1],
                                 protein_interaction_network.get_measurements(
                                     community, subset["time"],
                                     subset["post-translational modification"],
                                     combination.SITE_COMBINATION[subset.get(
                                         "site combination", "maxabs")])))

                            if configuration["community detection"][
                                    "Reactome enrichment"].get(
                                        "intersection", False):
                                cmty_proteins[community].intersection_update(
                                    protein_interaction_network.get_proteins(
                                        community,
                                        subset["time"],
                                        subset[
                                            "post-translational modification"],
                                        combination.SITE_COMBINATION[subset.get(
                                            "site combination", "maxabs")],
                                        lambda measurement, measurement_range=
                                        measurement_range: measurement <=
                                        measurement_range[0] or measurement >=
                                        measurement_range[1]))

                            else:
                                cmty_proteins[community].update(
                                    protein_interaction_network.get_proteins(
                                        community,
                                        subset["time"],
                                        subset[
                                            "post-translational modification"],
                                        combination.SITE_COMBINATION[subset.get(
                                            "site combination", "maxabs")],
                                        lambda measurement, measurement_range=
                                        measurement_range: measurement <=
                                        measurement_range[0] or measurement >=
                                        measurement_range[1]))

                reactome_enrichment = reactome.get_enrichment(
                    [
                        frozenset(cmty_proteins[community])
                        for community in communities
                    ],
                    enrichment_test=test.ENRICHMENT_TEST[
                        configuration["community detection"]
                        ["Reactome enrichment"].get("test", "hypergeometric")],
                    multiple_testing_correction=correction.CORRECTION[
                        configuration["community detection"]
                        ["Reactome enrichment"].get("correction",
                                                    "Benjamini-Hochberg")],
                    organism=configuration["community detection"]
                    ["Reactome enrichment"].get("organism", 9606),
                    annotation_as_reference=configuration["community detection"]
                    ["Reactome enrichment"].get("annotation", False))

                for k, community in enumerate(sorted(
                        communities,
                        key=lambda community: community.number_of_nodes(),
                        reverse=True),
                                              start=1):
                    for (pathway,
                         name), p in sorted(reactome_enrichment[frozenset(
                             cmty_proteins[community])].items(),
                                            key=lambda item: item[1]):
                        if p <= configuration["community detection"][
                                "Reactome enrichment"].get("p", 1.0):
                            export[community] = True
                            logger.info(f"{k}\t{pathway}\t{p:.2e}\t{name}")

            else:
                reactome_enrichment = reactome.get_enrichment(
                    [frozenset(community.nodes()) for community in communities],
                    enrichment_test=test.ENRICHMENT_TEST[
                        configuration["community detection"]
                        ["Reactome enrichment"].get("test", "hypergeometric")],
                    multiple_testing_correction=correction.CORRECTION[
                        configuration["community detection"]
                        ["Reactome enrichment"].get("correction",
                                                    "Benjamini-Hochberg")],
                    organism=configuration["community detection"]
                    ["Reactome enrichment"].get("organism", 9606),
                    annotation_as_reference=configuration["community detection"]
                    ["Reactome enrichment"].get("annotation", False))

                for k, community in enumerate(sorted(
                        communities,
                        key=lambda community: community.number_of_nodes(),
                        reverse=True),
                                              start=1):
                    for (pathway,
                         name), p in sorted(reactome_enrichment[frozenset(
                             community.nodes())].items(),
                                            key=lambda item: item[1]):
                        if p <= configuration["community detection"][
                                "Reactome enrichment"].get("p", 1.0):
                            export[community] = True
                            logger.info(f"{k}\t{pathway}\t{p:.2e}\t{name}")

        if "measurement enrichment" in configuration["community detection"]:
            enrichment = protein_interaction_network.get_measurement_enrichment(
                network,
                communities,
                measurements=default.MEASUREMENT_RANGE[
                    configuration["community detection"]
                    ["measurement enrichment"].get("conversion")],
                convert_measurement=conversion.MEASUREMENT_CONVERSION[
                    configuration["community detection"]
                    ["measurement enrichment"].get("conversion")],
                site_combination=combination.SITE_COMBINATION[
                    configuration["community detection"]
                    ["measurement enrichment"].get("site combination",
                                                   "maxabs")],
                enrichment_test=test.ENRICHMENT_TEST[
                    configuration["community detection"]
                    ["measurement enrichment"].get("test", "hypergeometric")],
                multiple_testing_correction=correction.CORRECTION[
                    configuration["community detection"]
                    ["measurement enrichment"].get("correction",
                                                   "Benjamini-Hochberg")])

            for k, community in enumerate(sorted(
                    communities,
                    key=lambda community: community.number_of_nodes(),
                    reverse=True),
                                          start=1):
                for time in enrichment[community]:
                    for modification, p in sorted(
                            enrichment[community][time].items(),
                            key=lambda item: item[1]):
                        if p <= configuration["community detection"][
                                "measurement enrichment"].get("p", 1.0):
                            export[community] = True
                            logger.info(f"{k}\tmeasurement enrichment\t{time} "
                                        f"{modification}\t{p:.2e}")

        if "measurement location" in configuration["community detection"]:
            location = protein_interaction_network.get_measurement_location(
                network,
                communities,
                combination.SITE_COMBINATION[
                    configuration["community detection"]
                    ["measurement location"].get("site combination", "maxabs")],
                location_test=test.LOCATION_TEST[
                    configuration["community detection"]
                    ["measurement location"].get("test", "Wilcoxon")],
                multiple_testing_correction=correction.CORRECTION[
                    configuration["community detection"]
                    ["measurement location"].get("correction",
                                                 "Benjamini-Hochberg")])

            for k, community in enumerate(sorted(
                    communities,
                    key=lambda community: community.number_of_nodes(),
                    reverse=True),
                                          start=1):
                for time in location[community]:
                    for modification, p in sorted(
                            location[community][time].items(),
                            key=lambda item: item[1]):
                        if p <= configuration["community detection"][
                                "measurement location"].get("p", 1.0):
                            export[community] = True
                            logger.info(f"{k}\tmeasurement location\t{time} "
                                        f"{modification}\t{p:.2e}")

        for k, community in enumerate(sorted(
                communities,
                key=lambda community: community.number_of_nodes(),
                reverse=True),
                                      start=1):
            if export[community]:
                protein_interaction_network.export(
                    community,
                    f"{logger.name}.{index}.{k}"
                    if index else f"{logger.name}.{k}",
                )


def process_configuration(configurations: Sequence[Mapping[str, Any]],
                          logger: logging.Logger) -> None:
    """
    Executes workflows specified in configurations sequentially.

    Args:
        configurations: The specification of workflows.
        logger: A configuration-specific logger.
    """
    if len(configurations) > 1:
        for i, configuration in enumerate(configurations, start=1):
            process_workflow(configuration, logger, i)
    else:
        process_workflow(configurations[0], logger)


def process_configuration_file(configuration_file: str) -> None:
    """
    Launches execution of a workflow.

    Args:
        configuration_file: file name of configuration file.
    """
    with open(configuration_file, encoding="utf-8") as configuration:
        process_configuration(
            json.load(configuration),
            logging.getLogger(
                os.path.splitext(os.path.basename(configuration_file))[0]))


def main() -> None:
    """Launches workflow execution."""
    parser = argparse.ArgumentParser()

    parser.add_argument("-c",
                        "--configurations",
                        help="concurrently processed configuration files each "
                        "specifying an array of workflows processed "
                        "sequentially",
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