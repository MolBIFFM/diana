"""
DIANA: data integration and network analysis of post-translational modification
based on mass spectrometry data
"""
import argparse
import concurrent.futures
import csv
import json
import logging
import os
import re
import sys
from typing import Any, Mapping

import networkx as nx

from cytoscape import (gene_ontology_network_style,
                       protein_interaction_network_style,
                       reactome_network_style)
from databases import gene_ontology, reactome
from interface import (average, score, correction, default, modularization,
                       order, prioritization, test)
from networks import (gene_ontology_network, protein_interaction_network,
                      reactome_network)


def process_workflow(identifier: str, configuration: Mapping[str, Any]) -> None:
    """
    Executes a workflow with identifier specified in configuration.

    Args:
        identifier: An identifier for the workflow.
        configuration: The specification of a workflow.
    """
    # Construct a workflow-specific logger.
    logger = logging.getLogger(identifier)

    # Initialize the protein-protein interaction network.
    network = protein_interaction_network.get_network()

    # Incorporate protein for mass spectrometry data sets into the
    # protein-protein interaction network.
    for time in configuration.get("mass spectrometry", {}):
        for modification in configuration["mass spectrometry"][time]:
            if not (configuration["mass spectrometry"][time][modification].get(
                    "file") and os.path.isfile(
                        configuration["mass spectrometry"][time][modification]
                        ["file"]) and configuration["mass spectrometry"][time]
                    [modification].get("accession column")):
                if not configuration["mass spectrometry"][time][
                        modification].get("file"):
                    logger.warning(
                        "File for modification %s at time %s is not specified.",
                        modification, time)

                elif not os.path.isfile(configuration["mass spectrometry"][time]
                                        [modification]["file"]):
                    logger.warning(
                        "File specified for modification %s at time %s does "
                        "not exist.", modification, time)

                elif not configuration["mass spectrometry"][time][
                        modification].get("accession column"):
                    logger.warning(
                        "Accession column for modification %s at time %s is "
                        "not specified.", modification, time)

                continue

            logger.info(
                "Adding proteins for modification %s at time %s from %s.",
                modification, time,
                configuration["mass spectrometry"][time][modification]["file"])

            # Incorporate proteins associated with site-specific measurements
            # into the protein-protein interaction network.
            if configuration["mass spectrometry"][time][modification].get(
                    "position column"):
                protein_interaction_network.add_sites_from_table(
                    network,
                    file=configuration["mass spectrometry"][time][modification]
                    ["file"],
                    protein_accession_column=configuration["mass spectrometry"]
                    [time][modification]["accession column"],
                    protein_accession_format=re.compile(
                        configuration["mass spectrometry"][time]
                        [modification].get("accession format", "^(.+)$")),
                    position_column=configuration["mass spectrometry"][time]
                    [modification]["position column"],
                    position_format=re.compile(
                        configuration["mass spectrometry"][time]
                        [modification].get("position format", "^(.+)$")),
                    replicate_columns=configuration["mass spectrometry"][time]
                    [modification]["replicate columns"],
                    replicate_format=re.compile(
                        configuration["mass spectrometry"][time]
                        [modification].get("replicate format", "^(.+)$")),
                    sheet_name=configuration["mass spectrometry"][time]
                    [modification].get("sheet", 1) - 1 if isinstance(
                        configuration["mass spectrometry"][time]
                        [modification].get("sheet", 1),
                        int) else configuration["mass spectrometry"][time]
                    [modification]["sheet"],
                    header=configuration["mass spectrometry"][time]
                    [modification].get("header", 1) - 1,
                    time=int(time) if time.isnumeric() else 0,
                    modification=modification,
                    number_sites=configuration["mass spectrometry"][time]
                    [modification].get("sites", 5),
                    number_replicates=configuration["mass spectrometry"][time]
                    [modification].get("replicates", 1),
                    replicate_average=average.REPLICATE_AVERAGE[
                        configuration["mass spectrometry"][time]
                        [modification].get("replicate average", "mean")],
                    measurement_score=score.LOGARITHM[
                        configuration["mass spectrometry"][time]
                        [modification].get("logarithm")],
                    site_prioritization=prioritization.SITE_PRIORITIZATION[
                        configuration["mass spectrometry"][time]
                        [modification].get("site prioritization", "absolute")],
                    site_order=order.SITE_ORDER[
                        configuration["mass spectrometry"][time]
                        [modification].get("site order", "measurement")])

            # Incorporate proteins associated with protein-specific
            # measurements into the protein-protein interaction network.
            else:
                protein_interaction_network.add_proteins_from_table(
                    network,
                    file=configuration["mass spectrometry"][time][modification]
                    ["file"],
                    protein_accession_column=configuration["mass spectrometry"]
                    [time][modification]["accession column"],
                    protein_accession_format=re.compile(
                        configuration["mass spectrometry"][time]
                        [modification].get("accession format", "^(.+)$")),
                    replicate_columns=configuration["mass spectrometry"][time]
                    [modification]["replicate columns"],
                    replicate_format=re.compile(
                        configuration["mass spectrometry"][time]
                        [modification].get("replicate format", "^(.+)$")),
                    sheet_name=configuration["mass spectrometry"][time]
                    [modification].get("sheet", 1) - 1 if isinstance(
                        configuration["mass spectrometry"][time]
                        [modification].get("sheet", 1),
                        int) else configuration["mass spectrometry"][time]
                    [modification]["sheet"],
                    header=configuration["mass spectrometry"][time]
                    [modification].get("header", 1) - 1,
                    time=int(time) if time.isnumeric() else 0,
                    modification=modification,
                    number_replicates=configuration["mass spectrometry"][time]
                    [modification].get("replicates", 1),
                    measurement_score=score.LOGARITHM[
                        configuration["mass spectrometry"][time]
                        [modification].get("logarithm")])

    # Incorporate protein-protein interaction networks into the protein-protein
    # interaction network.
    for i, item in enumerate(configuration.get("networks", {})):
        if not (item.get("network") and os.path.isfile(item["network"])):
            if not item.get("network"):
                logger.warning(
                    "File for protein-protein interaction network %d is not "
                    "specified.", i)

            elif not os.path.isfile(item["network"]):
                logger.warning(
                    "File for protein-protein interaction network %d does not "
                    "exist.", i)

            continue

        logger.info("Adding proteins from %s.", item["network"])

        network = nx.compose(network,
                             nx.readwrite.graphml.read_graphml(item["network"]))

    # Incorporate protein accessions into the protein-protein interaction
    # network.
    for item in configuration.get("proteins", {}):
        if item.get("accessions"):
            logger.info("Adding proteins from accessions.")

            network.add_nodes_from(
                protein_accession for protein_accession in item["accessions"]
                # https://www.uniprot.org/help/accession_numbers
                # https://www.uniprot.org/help/alternative_products
                if re.fullmatch(
                    r"([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]"
                    r"([A-Z][A-Z0-9]{2}[0-9]){1,2})(-[1-9][0-9]+)?",
                    protein_accession))

    # Assess availability of a local file for UniProt protein accessions.
    if configuration.get("UniProt", {}).get("file"):
        if os.path.isfile(configuration["UniProt"]["file"]):
            logger.info("Parsing UniProt protein accessions from %s.",
                        configuration["UniProt"]["file"])
        else:
            logger.warning(
                "File specified for UniProt protein accessions does not exist.")

    # Map protein accessions to primary UniProt accessions and remove proteins
    # not found in Swiss-Prot from the protein-protein interaction network.
    nodes_to_remove = set(network.nodes())
    for organism in set(
            configuration.get("mass spectrometry", {}).get(time, {}).get(
                modification, {}).get("organism", 9606)
            for time in configuration.get("mass spectrometry", {})
            for modification in configuration.get("mass spectrometry", {}).get(
                time, {})) | set(
                    item.get("organism", 9606)
                    for item in configuration.get("networks", {})) | set(
                        item.get("organism", 9606)
                        for item in configuration.get("proteins", {})):
        nodes_to_remove.intersection_update(
            protein_interaction_network.map_proteins(network,
                                                     organism,
                                                     file=configuration.get(
                                                         "UniProt",
                                                         {}).get("file")))
    network.remove_nodes_from(nodes_to_remove)

    # Add neighboring proteins to the protein-protein interaction network.
    if "protein-protein interactions" in configuration:
        for neighbors in range(
                max(configuration["protein-protein interactions"].get(
                    database, {}).get("neighbors", 0) for database in
                    configuration["protein-protein interactions"])):
            interacting_proteins = set()

            # Add neighboring proteins from BioGRID to the protein-protein
            # interaction network.
            if "BioGRID" in configuration[
                    "protein-protein interactions"] and configuration[
                        "protein-protein interactions"]["BioGRID"].get(
                            "neighbors", 0) > neighbors:
                if configuration["protein-protein interactions"]["BioGRID"].get(
                        "file") and os.path.isfile(
                            configuration["protein-protein interactions"]
                            ["BioGRID"]["file"]):
                    logger.info(
                        "Adding neighbors of order %d from BioGRID from %s.",
                        neighbors, configuration["protein-protein interactions"]
                        ["BioGRID"]["file"])
                else:
                    if configuration["protein-protein interactions"][
                            "BioGRID"].get("file"):
                        logger.warning(
                            "File specified for protein-protein interactions "
                            "from BioGRID does not exist.")

                    logger.info("Adding neighbors of order %d from BioGRID.",
                                neighbors)

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
                        ["BioGRID"].get("version"),
                        file=configuration["protein-protein interactions"]
                        ["BioGRID"].get("file"),
                        file_uniprot=configuration.get("UniProt",
                                                       {}).get("file")))

            # Add neighboring proteins from CORUM to the protein-protein
            # interaction network.
            if "CORUM" in configuration[
                    "protein-protein interactions"] and configuration[
                        "protein-protein interactions"]["CORUM"].get(
                            "neighbors", 0) > neighbors:
                if configuration["protein-protein interactions"]["CORUM"].get(
                        "file") and os.path.isfile(
                            configuration["protein-protein interactions"]
                            ["CORUM"]["file"]):
                    logger.info(
                        "Adding neighbors of order %d from CORUM from %s.",
                        neighbors, configuration["protein-protein interactions"]
                        ["CORUM"]["file"])
                else:
                    if configuration["protein-protein interactions"][
                            "CORUM"].get("file"):
                        logger.warning(
                            "File specified for protein-protein interactions "
                            "from CORUM does not exist.")

                    logger.info("Adding neighbors of order %d from CORUM.",
                                neighbors)

                interacting_proteins.update(
                    protein_interaction_network.get_neighbors_from_corum(
                        network,
                        purification_methods=configuration[
                            "protein-protein interactions"]["CORUM"].get(
                                "purification methods", []),
                        file=configuration["protein-protein interactions"]
                        ["CORUM"].get("file"),
                        file_uniprot=configuration.get("UniProt",
                                                       {}).get("file")))

            # Add neighboring proteins from IntAct to the protein-protein
            # interaction network.
            if "IntAct" in configuration[
                    "protein-protein interactions"] and configuration[
                        "protein-protein interactions"]["IntAct"].get(
                            "neighbors", 0) > neighbors:
                if configuration["protein-protein interactions"]["IntAct"].get(
                        "file") and os.path.isfile(
                            configuration["protein-protein interactions"]
                            ["IntAct"]["file"]):
                    logger.info(
                        "Adding neighbors of order %d from IntAct from %s.",
                        neighbors, configuration["protein-protein interactions"]
                        ["IntAct"]["file"])
                else:
                    if configuration["protein-protein interactions"][
                            "IntAct"].get("file"):
                        logger.warning(
                            "File specified for protein-protein interactions "
                            "from IntAct does not exist.")

                    logger.info("Adding neighbors of order %d from IntAct.",
                                neighbors)

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
                        file=configuration["protein-protein interactions"]
                        ["IntAct"].get("file"),
                        file_uniprot=configuration.get("UniProt",
                                                       {}).get("file")))

            # Add neighboring proteins from MINT to the protein-protein
            # interaction network.
            if "MINT" in configuration[
                    "protein-protein interactions"] and configuration[
                        "protein-protein interactions"]["MINT"].get(
                            "neighbors", 0) > neighbors:
                if configuration["protein-protein interactions"]["MINT"].get(
                        "file") and os.path.isfile(
                            configuration["protein-protein interactions"]
                            ["MINT"]["file"]):
                    logger.info(
                        "Adding neighbors of order %d from MINT from %s.",
                        neighbors, configuration["protein-protein interactions"]
                        ["MINT"]["file"])
                else:
                    if configuration["protein-protein interactions"][
                            "MINT"].get("file"):
                        logger.warning(
                            "File specified for protein-protein interactions "
                            "from MINT does not exist.")

                    logger.info("Adding neighbors of order %d from MINT.",
                                neighbors)

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
                        file=configuration["protein-protein interactions"]
                        ["MINT"].get("file"),
                        file_uniprot=configuration.get("UniProt",
                                                       {}).get("file")))

            # Add neighboring proteins from Reactome to the protein-protein
            # interaction network.
            if "Reactome" in configuration[
                    "protein-protein interactions"] and configuration[
                        "protein-protein interactions"]["Reactome"].get(
                            "neighbors", 0) > neighbors:
                if configuration["protein-protein interactions"][
                        "Reactome"].get("file") and os.path.isfile(
                            configuration["protein-protein interactions"]
                            ["Reactome"]["file"]):
                    logger.info(
                        "Adding neighbors of order %d from Reactome from %s.",
                        neighbors, configuration["protein-protein interactions"]
                        ["Reactome"]["file"])
                else:
                    if configuration["protein-protein interactions"][
                            "Reactome"].get("file"):
                        logger.warning(
                            "File specified for protein-protein interactions "
                            "from Reactome does not exist.")

                    logger.info("Adding neighbors of order %d from Reactome.",
                                neighbors)

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
                        file=configuration["protein-protein interactions"]
                        ["Reactome"].get("file"),
                        file_uniprot=configuration.get("UniProt",
                                                       {}).get("file")))

            # Add neighboring proteins from STRING to the protein-protein
            # interaction network.
            if "STRING" in configuration[
                    "protein-protein interactions"] and configuration[
                        "protein-protein interactions"]["STRING"].get(
                            "neighbors", 0) > neighbors:
                if configuration["protein-protein interactions"]["STRING"].get(
                        "file") and os.path.isfile(
                            configuration["protein-protein interactions"]
                            ["STRING"]["file"]):
                    logger.info(
                        "Adding neighbors of order %d from STRING from %s.",
                        neighbors, configuration["protein-protein interactions"]
                        ["STRING"]["file"])
                else:
                    if configuration["protein-protein interactions"][
                            "STRING"].get("file"):
                        logger.warning(
                            "File specified for protein-protein interactions "
                            "from STRING does not exist.")

                    logger.info("Adding neighbors of order %d from STRING.",
                                neighbors)

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
                        any_score=configuration["protein-protein interactions"]
                        ["STRING"].get("any score", False),
                        file=configuration["protein-protein interactions"]
                        ["STRING"].get("file"),
                        file_accession_map=configuration[
                            "protein-protein interactions"]["STRING"].get(
                                "file accession map"),
                        file_uniprot=configuration.get("UniProt",
                                                       {}).get("file")))

            # Map protein accessions of neighboring proteins to primary UniProt
            # accessions and remove neighboring proteins not found in Swiss-Prot
            # from the protein-protein interaction network.
            network.add_nodes_from(interacting_proteins)
            nodes_to_remove = set(network.nodes())
            for organism in set(configuration["protein-protein interactions"]
                                [database].get("organism", 9606) for database in
                                configuration["protein-protein interactions"]):
                nodes_to_remove.intersection_update(
                    protein_interaction_network.map_proteins(
                        network, organism,
                        configuration.get("UniProt", {}).get("file")))
            network.remove_nodes_from(nodes_to_remove)

        # Add protein-protein interactions from BioGRID to the protein-protein
        # interaction network.
        if "BioGRID" in configuration["protein-protein interactions"]:
            if configuration["protein-protein interactions"]["BioGRID"].get(
                    "file") and os.path.isfile(
                        configuration["protein-protein interactions"]["BioGRID"]
                        ["file"]):
                logger.info(
                    "Adding protein-protein interactions from BioGRID from %s.",
                    configuration["protein-protein interactions"]["BioGRID"]
                    ["file"])
            else:
                if configuration["protein-protein interactions"]["BioGRID"].get(
                        "file"):
                    logger.warning(
                        "File specified for protein-protein interactions from "
                        "BioGRID does not exist.")

                logger.info("Adding protein-protein interactions from BioGRID.")

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
                ["BioGRID"].get("version"),
                file=configuration["protein-protein interactions"]
                ["BioGRID"].get("file"),
                file_uniprot=configuration.get("UniProt", {}).get("file"))

        # Add protein-protein interactions from CORUM to the protein-protein
        # interaction network.
        if "CORUM" in configuration["protein-protein interactions"]:
            if configuration["protein-protein interactions"]["CORUM"].get(
                    "file") and os.path.isfile(
                        configuration["protein-protein interactions"]["CORUM"]
                        ["file"]):
                logger.info(
                    "Adding protein-protein interactions from CORUM from %s.",
                    configuration["protein-protein interactions"]["CORUM"]
                    ["file"])
            else:
                if configuration["protein-protein interactions"]["CORUM"].get(
                        "file"):
                    logger.warning(
                        "File specified for protein-protein interactions from "
                        "CORUM does not exist.")

                logger.info("Adding protein-protein interactions from CORUM.")

            protein_interaction_network.add_protein_interactions_from_corum(
                network,
                purification_methods=configuration[
                    "protein-protein interactions"]["CORUM"].get(
                        "purification methods", []),
                file=configuration["protein-protein interactions"]["CORUM"].get(
                    "file"),
                file_uniprot=configuration.get("UniProt", {}).get("file"))

        # Add protein-protein interactions from IntAct to the protein-protein
        # interaction network.
        if "IntAct" in configuration["protein-protein interactions"]:
            if configuration["protein-protein interactions"]["IntAct"].get(
                    "file") and os.path.isfile(
                        configuration["protein-protein interactions"]["IntAct"]
                        ["file"]):
                logger.info(
                    "Adding protein-protein interactions from IntAct from %s.",
                    configuration["protein-protein interactions"]["IntAct"]
                    ["file"])
            else:
                if configuration["protein-protein interactions"]["IntAct"].get(
                        "file"):
                    logger.warning(
                        "File specified for protein-protein interactions from "
                        "IntAct does not exist.")

                logger.info("Adding protein-protein interactions from IntAct.")

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
                file=configuration["protein-protein interactions"]
                ["IntAct"].get("file"),
                file_uniprot=configuration.get("UniProt", {}).get("file"))

        # Add protein-protein interactions from MINT to the protein-protein
        # interaction network.
        if "MINT" in configuration["protein-protein interactions"]:
            if configuration["protein-protein interactions"]["MINT"].get(
                    "file") and os.path.isfile(
                        configuration["protein-protein interactions"]["MINT"]
                        ["file"]):
                logger.info(
                    "Adding protein-protein interactions from MINT from %s.",
                    configuration["protein-protein interactions"]["MINT"]
                    ["file"])
            else:
                if configuration["protein-protein interactions"]["MINT"].get(
                        "file"):
                    logger.warning(
                        "File specified for protein-protein interactions from "
                        "MINT does not exist.")

                logger.info("Adding protein-protein interactions from MINT.")

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
                file=configuration["protein-protein interactions"]["MINT"].get(
                    "file"),
                file_uniprot=configuration.get("UniProt", {}).get("file"))

        # Add protein-protein interactions from Reactome to the protein-protein
        # interaction network.
        if "Reactome" in configuration["protein-protein interactions"]:
            if configuration["protein-protein interactions"]["Reactome"].get(
                    "file") and os.path.isfile(
                        configuration["protein-protein interactions"]
                        ["Reactome"]["file"]):
                logger.info(
                    "Adding protein-protein interactions from Reactome from "
                    "%s.", configuration["protein-protein interactions"]
                    ["Reactome"]["file"])
            else:
                if configuration["protein-protein interactions"][
                        "Reactome"].get("file"):
                    logger.warning(
                        "File specified for protein-protein interactions from "
                        "Reactome does not exist.")

                logger.info(
                    "Adding protein-protein interactions from Reactome.")

            protein_interaction_network.add_protein_interactions_from_reactome(
                network,
                interaction_context=configuration[
                    "protein-protein interactions"]["Reactome"].get(
                        "interaction context", []),
                interaction_type=configuration["protein-protein interactions"]
                ["Reactome"].get("interaction type", []),
                organism=configuration["protein-protein interactions"]
                ["Reactome"].get("organism", 9606),
                file=configuration["protein-protein interactions"]
                ["Reactome"].get("file"),
                file_uniprot=configuration.get("UniProt", {}).get("file"))

        # Add protein-protein interactions from STRING to the protein-protein
        # interaction network.
        if "STRING" in configuration["protein-protein interactions"]:
            if configuration["protein-protein interactions"]["STRING"].get(
                    "file") and os.path.isfile(
                        configuration["protein-protein interactions"]["STRING"]
                        ["file"]):
                logger.info(
                    "Adding protein-protein interactions from STRING from %s.",
                    configuration["protein-protein interactions"]["STRING"]
                    ["file"])
            else:
                if configuration["protein-protein interactions"]["STRING"].get(
                        "file"):
                    logger.warning(
                        "File specified for protein-protein interactions from "
                        "STRING does not exist.")

                logger.info("Adding protein-protein interactions from STRING.")

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
                any_score=configuration["protein-protein interactions"]
                ["STRING"].get("any score", False),
                file=configuration["protein-protein interactions"]
                ["STRING"].get("file"),
                file_accession_map=configuration["protein-protein interactions"]
                ["STRING"].get("file accession map"),
                file_uniprot=configuration.get("UniProt", {}).get("file"))

        # Compile Cytoscape styles for the protein-protein interaction network
        # and annotate the protein-protein interaction network.
        if any(
                protein_interaction_network.get_modifications(network, time)
                for time in protein_interaction_network.get_times(network)):

            if "Cytoscape" in configuration:
                # Annotate nodes of the protein-protein interaction network with
                # a categorization of average measurements for different types
                # of post-translational modification at different times of
                # measurement.
                measurements = {
                    modification: default.MEASUREMENT_RANGE.get(
                        score, default.MEASUREMENT_RANGE[None])
                    for modification, score in configuration["Cytoscape"].get(
                        "node color", {}).get("score", {}).items()
                }

                for modification, measurement_range in configuration[
                        "Cytoscape"].get("node color",
                                         {}).get("measurement", {}).items():
                    measurements[modification] = measurement_range

                protein_interaction_network.set_measurements(
                    network,
                    site_average={
                        modification: average.SITE_AVERAGE.get(
                            site_average,
                            average.SITE_AVERAGE["maximum absolute logarithm"])
                        for modification,
                        site_average in configuration["Cytoscape"].get(
                            "site average", {}).items()
                    },
                    replicate_average={
                        modification: average.REPLICATE_AVERAGE.get(
                            replicate_average,
                            average.REPLICATE_AVERAGE["mean"]) for modification,
                        replicate_average in configuration["Cytoscape"].get(
                            "replicate average", {}).items()
                    },
                    measurements=measurements,
                    measurement_score={
                        modification: score.MEASUREMENT_SCORE.get(
                            measurement_score, score.MEASUREMENT_SCORE[None])
                        for modification, measurement_score in
                        configuration["Cytoscape"].get("score", {})
                    })

                # Annotate edges of the protein-protein interaction network with
                # average protein-protein interaction confidence score reflected
                # as edge transparency in Cytoscape.
                protein_interaction_network.set_edge_weights(
                    network,
                    weight=average.CONFIDENCE_SCORE_AVERAGE[
                        configuration["Cytoscape"].get("edge transparency")],
                    attribute="score")

                # Compile Cytoscape styles for the protein-protein interaction
                # network.
                styles = protein_interaction_network_style.get_styles(
                    network,
                    node_shape_modifications=configuration["Cytoscape"].get(
                        "node shape",
                        {}).get("post-translational modifications", []),
                    node_color_modifications=configuration["Cytoscape"].get(
                        "node color",
                        {}).get("post-translational modifications", []),
                    node_size_modification=configuration["Cytoscape"].get(
                        "node size", {}).get("post-translational modification",
                                             None),
                    bar_chart_modifications=configuration["Cytoscape"].get(
                        "bar chart", {}).get("post-translational modifications",
                                             []),
                    measurement_score={
                        modification: score.MEASUREMENT_SCORE.get(
                            measurement_score, score.MEASUREMENT_SCORE[None])
                        for modification, measurement_score in
                        configuration["Cytoscape"].get("score", {})
                    },
                    site_average={
                        modification: average.SITE_AVERAGE.get(
                            site_average,
                            average.SITE_AVERAGE["maximum absolute logarithm"])
                        for modification,
                        site_average in configuration["Cytoscape"].get(
                            "site average", {}).items()
                    },
                    replicate_average={
                        modification: average.REPLICATE_AVERAGE.get(
                            replicate_average,
                            average.REPLICATE_AVERAGE["mean"]) for modification,
                        replicate_average in configuration["Cytoscape"].get(
                            "replicate average", {}).items()
                    },
                    confidence_score_average=average.CONFIDENCE_SCORE_AVERAGE[
                        configuration["Cytoscape"].get("edge transparency")])

                # Export the Cytoscape styles.
                file = protein_interaction_network_style.export(
                    styles, identifier)

                if file is None:
                    logger.warning(
                        "The Cytoscape styles for the protein-protein "
                        "interaction network were not exported due to naming "
                        "conflict.")
                else:
                    logger.info(
                        "The Cytoscape styles for the protein-protein "
                        "interaction network were exported to %s.", file)

        # Export the protein-protein interaction network.
        file = protein_interaction_network.export(network, identifier)

        if file is None:
            logger.warning(
                "The protein-protein interaction network was not exported due "
                "to naming conflict.")
        else:
            logger.info(
                "The protein-protein interaction network was exported to %s.",
                file)

    # Assess availability of a local file for the Gene Ontology.
    if configuration.get("Gene Ontology", {}).get("file", {}).get("ontology"):
        if os.path.isfile(configuration["Gene Ontology"]["file"]["ontology"]):
            logger.info("Parsing Gene Ontology from %s.",
                        configuration["Gene Ontology"]["file"]["ontology"])
        else:
            logger.warning("File specified for Gene Ontology does not exist.")

    # Assess availability of a local file for Gene Ontology protein annotations.
    if configuration.get("Gene Ontology", {}).get("file", {}).get("annotation"):
        if os.path.isfile(configuration["Gene Ontology"]["file"]["annotation"]):
            logger.info("Parsing Gene Ontology annotation from %s.",
                        configuration["Gene Ontology"]["file"]["annotation"])
        else:
            logger.warning(
                "File specified for Gene Ontology annotation does not exist.")

    # Assess availability of a local file for Gene Ontology protein isoform
    # annotations.
    if configuration.get("Gene Ontology", {}).get("file",
                                                  {}).get("annotation isoform"):
        if os.path.isfile(
                configuration["Gene Ontology"]["file"]["annotation isoform"]):
            logger.info(
                "Parsing Gene Ontology isoform annotation from %s.",
                configuration["Gene Ontology"]["file"]["annotation isoform"])
        else:
            logger.warning(
                "File specified for Gene Ontology isoform annotation does not "
                "exist.")

    # Assess availability of a local file for Reactome pathways.
    if configuration.get("Reactome", {}).get("file", {}).get("pathways"):
        if os.path.isfile(configuration["Reactome"]["file"]["pathways"]):
            logger.info("Parsing Reactome pathways from %s.",
                        configuration["Reactome"]["file"]["pathways"])
        else:
            logger.warning(
                "File specified for Reactome pathways does not exist.")

    # Assess availability of a local file for Reactome pathway relations.
    if configuration.get("Reactome", {}).get("file",
                                             {}).get("pathways relation"):
        if os.path.isfile(
                configuration["Reactome"]["file"]["pathways relation"]):
            logger.info("Parsing Reactome pathway relations from %s.",
                        configuration["Reactome"]["file"]["pathways relation"])
        else:
            logger.warning(
                "File specified for Reactome pathway relations does not exist.")

    # Assess availability of a local file for associations of Reactome pathways
    # and UniProt protein accessions.
    if configuration.get("Reactome", {}).get("file", {}).get("accession map"):
        if os.path.isfile(configuration["Reactome"]["file"]["accession map"]):
            logger.info(
                "Parsing associations of Reactome pathways and UniProt protein "
                "accessions from %s.",
                configuration["Reactome"]["file"]["accession map"])
        else:
            logger.warning(
                "File specified associations of Reactome pathways and UniProt "
                "protein accessions does not exist.")

    # Export a Gene Ontology network.
    if "Gene Ontology network" in configuration:
        logger.info("Compiling the Gene Ontology network.")

        # Compile a Gene Ontology network from subsets of proteins represented
        # in the protein-protein interaction network determined from measurement
        # averages.
        if "post-translational modifications" in configuration[
                "Gene Ontology network"]:
            if configuration["Gene Ontology network"].get(
                    "intersection", False):
                proteins = set(network.nodes())
            else:
                proteins = set()

            for time in protein_interaction_network.get_times(network):
                for m in configuration["Gene Ontology network"].get(
                        "post-translational modifications", []):
                    if m in protein_interaction_network.get_modifications(
                            network, time):
                        measurement_range = (
                            score.MEASUREMENT_SCORE[
                                configuration["Gene Ontology network"].get(
                                    "score", {}).get(m)]
                            (configuration["Gene Ontology network"].
                             get("measurement", {}).get(
                                 m, default.MEASUREMENT_RANGE[
                                     configuration["Gene Ontology network"].get(
                                         "score", {}).get(m)])[0],
                             protein_interaction_network.get_measurements(
                                 network, time, m, average.SITE_AVERAGE[
                                     configuration["Gene Ontology network"].get(
                                         "site average",
                                         {}).get(m,
                                                 "maximum absolute logarithm")],
                                 average.REPLICATE_AVERAGE[
                                     configuration["Gene Ontology network"].get(
                                         "replicate average",
                                         {}).get(m, "mean")])),
                            score.MEASUREMENT_SCORE[
                                configuration["Gene Ontology network"].get(
                                    "score", {}).get(m)]
                            (configuration["Gene Ontology network"].get(
                                "measurement", {}).get(
                                    m, default.MEASUREMENT_RANGE[configuration[
                                        "Gene Ontology network"].get(
                                            "score", {}).get(m)])[1],
                             protein_interaction_network.get_measurements(
                                 network, time, m, average.SITE_AVERAGE[
                                     configuration["Gene Ontology network"].get(
                                         "site average",
                                         {}).get(m,
                                                 "maximum absolute logarithm")],
                                 average.REPLICATE_AVERAGE[
                                     configuration["Gene Ontology network"].get(
                                         "replicate average",
                                         {}).get(m, "mean")])))

                        if configuration["Gene Ontology network"].get(
                                "intersection", False):
                            proteins.intersection_update(
                                protein_interaction_network.get_proteins(
                                    network, time, m,
                                    average.SITE_AVERAGE[configuration[
                                        "Gene Ontology network"].get(
                                            "site average", {}).get(
                                                m,
                                                "maximum absolute logarithm")],
                                    average.REPLICATE_AVERAGE[configuration[
                                        "Gene Ontology network"].get(
                                            "replicate average",
                                            {}).get(m, "mean")],
                                    measurement_range))

                        else:
                            proteins.update(
                                protein_interaction_network.get_proteins(
                                    network, time, m,
                                    average.SITE_AVERAGE[configuration[
                                        "Gene Ontology network"].get(
                                            "site average", {}).get(
                                                m,
                                                "maximum absolute logarithm")],
                                    average.REPLICATE_AVERAGE[configuration[
                                        "Gene Ontology network"].get(
                                            "replicate average",
                                            {}).get(m, "mean")],
                                    measurement_range))

            ontology_network = gene_ontology_network.get_network(
                proteins,
                reference=network.nodes()
                if not configuration["Gene Ontology network"].get(
                    "annotation", False) else None,
                namespaces=[
                    namespace.replace(" ", "_")
                    for namespace in configuration["Gene Ontology network"].get(
                        "namespaces", [])
                ],
                enrichment_test=test.ENRICHMENT_TEST[(
                    configuration["Gene Ontology network"].get(
                        "test", "hypergeometric"),
                    configuration["Gene Ontology network"].get(
                        "increase", True))],
                multiple_testing_correction=correction.CORRECTION[
                    configuration["Gene Ontology network"].get(
                        "correction", "Benjamini-Yekutieli")],
                organism=configuration["Gene Ontology network"].get(
                    "organism", 9606),
                file_ontology=configuration.get("Gene Ontology",
                                                {}).get("file",
                                                        {}).get("ontology"),
                file_annotation=configuration.get("Gene Ontology",
                                                  {}).get("file",
                                                          {}).get("annotation"),
                file_annotation_isoform=configuration.get(
                    "Gene Ontology", {}).get("file",
                                             {}).get("annotation isoform"),
                file_uniprot=configuration.get("UniProt", {}).get("file"))

        # Compile a Gene Ontology network from proteins represented in the
        # protein-protein interaction network.
        else:
            ontology_network = gene_ontology_network.get_network(
                network.nodes(),
                namespaces=[
                    namespace.replace(" ", "_")
                    for namespace in configuration["Gene Ontology network"].get(
                        "namespaces", [])
                ],
                enrichment_test=test.ENRICHMENT_TEST[(
                    configuration["Gene Ontology network"].get(
                        "test", "hypergeometric"),
                    configuration["Gene Ontology network"].get(
                        "increase", True))],
                multiple_testing_correction=correction.CORRECTION[
                    configuration["Gene Ontology network"].get(
                        "correction", "Benjamini-Yekutieli")],
                organism=configuration["Gene Ontology network"].get(
                    "organism", 9606),
                file_ontology=configuration.get("Gene Ontology",
                                                {}).get("file",
                                                        {}).get("ontology"),
                file_annotation=configuration.get("Gene Ontology",
                                                  {}).get("file",
                                                          {}).get("annotation"),
                file_annotation_isoform=configuration.get(
                    "Gene Ontology", {}).get("file",
                                             {}).get("annotation isoform"),
                file_uniprot=configuration.get("UniProt", {}).get("file"))

        # Export the Gene Ontology network.
        file = gene_ontology_network.export(ontology_network,
                                            f"{identifier}_gene_ontology")

        if file is None:
            logger.warning(
                "The Gene Ontology network was not exported due to naming "
                "conflict.")
        else:
            logger.info("The Gene Ontology network was exported to %s.", file)

        # Export the Cytoscape style for the Gene Ontology network.
        if "Cytoscape" in configuration:
            ontology_network_styles = gene_ontology_network_style.get_styles(
                ontology_network)
            file = gene_ontology_network_style.export(
                ontology_network_styles, f"{identifier}_gene_ontology")

            if file is None:
                logger.warning(
                    "The Cytoscape style for the Gene Ontology network was not "
                    "exported due to naming conflict.")
            else:
                logger.info(
                    "The Cytoscape style for the Gene Ontology network was "
                    "exported to %s.", file)

    # Export a Reactome network.
    if "Reactome network" in configuration:
        logger.info("Compiling the Reactome network.")

        # Compile a Reactome network from subsets of proteins represented in the
        # protein-protein interaction network determined from measurement
        # averages.
        if "post-translational modifications" in configuration[
                "Reactome network"]:
            if configuration["Reactome network"].get("intersection", False):
                proteins = set(network.nodes())
            else:
                proteins = set()

            for time in protein_interaction_network.get_times(network):
                for m in configuration["Reactome network"].get(
                        "post-translational modifications", []):
                    if m in protein_interaction_network.get_modifications(
                            network, time):
                        measurement_range = (
                            score.MEASUREMENT_SCORE[
                                configuration["Reactome network"].get(
                                    "score", {}).get(m)]
                            (configuration["Reactome network"].get(
                                "measurement", {}).get(
                                    m, default.MEASUREMENT_RANGE[
                                        configuration["Reactome network"].get(
                                            "score", {}).get(m)])[0],
                             protein_interaction_network.get_measurements(
                                 network, time, m, average.SITE_AVERAGE[
                                     configuration["Reactome network"].get(
                                         "site average",
                                         {}).get(m,
                                                 "maximum absolute logarithm")],
                                 average.REPLICATE_AVERAGE[
                                     configuration["Reactome network"].get(
                                         "replicate average",
                                         {}).get(m, "mean")])),
                            score.MEASUREMENT_SCORE[
                                configuration["Reactome network"].get(
                                    "score", {}).get(m)]
                            (configuration["Reactome network"].get(
                                "measurement", {}).get(
                                    m, default.MEASUREMENT_RANGE[
                                        configuration["Reactome network"].get(
                                            "score", {}).get(m)])[1],
                             protein_interaction_network.get_measurements(
                                 network, time, m, average.SITE_AVERAGE[
                                     configuration["Reactome network"].get(
                                         "site average",
                                         {}).get(m,
                                                 "maximum absolute logarithm")],
                                 average.REPLICATE_AVERAGE[
                                     configuration["Reactome network"].get(
                                         "replicate average",
                                         {}).get(m, "mean")])))

                        if configuration["Reactome network"].get(
                                "intersection", False):
                            proteins.intersection_update(
                                protein_interaction_network.get_proteins(
                                    network, time, m, average.SITE_AVERAGE[
                                        configuration["Reactome network"].get(
                                            "site average", {}).get(
                                                m,
                                                "maximum absolute logarithm")],
                                    average.REPLICATE_AVERAGE[
                                        configuration["Reactome network"].get(
                                            "replicate average",
                                            {}).get(m, "mean")],
                                    measurement_range))

                        else:
                            proteins.update(
                                protein_interaction_network.get_proteins(
                                    network, time, m, average.SITE_AVERAGE[
                                        configuration["Reactome network"].get(
                                            "site average", {}).get(
                                                m,
                                                "maximum absolute logarithm")],
                                    average.REPLICATE_AVERAGE[
                                        configuration["Reactome network"].get(
                                            "replicate average",
                                            {}).get(m, "mean")],
                                    measurement_range))

            pathway_network = reactome_network.get_network(
                proteins,
                reference=network.nodes() if
                not configuration["Reactome network"].get("annotation", False)
                else None,
                enrichment_test=test.ENRICHMENT_TEST[(
                    configuration["Reactome network"].get(
                        "test", "hypergeometric"),
                    configuration["Reactome network"].get("increase", True))],
                multiple_testing_correction=correction.CORRECTION[
                    configuration["Reactome network"].get(
                        "correction", "Benjamini-Yekutieli")],
                organism=configuration["Reactome network"].get(
                    "organism", 9606),
                file_pathways=configuration.get("Reactome",
                                                {}).get("file",
                                                        {}).get("pathways"),
                file_pathways_relation=configuration.get("Reactome", {}).get(
                    "file", {}).get("pathways relation"),
                file_accession_map=configuration.get("Reactome", {}).get(
                    "file", {}).get("accession map"),
                file_uniprot=configuration.get("UniProt", {}).get("file"))

        # Compile a Reactome network from proteins represented in the
        # protein-protein interaction network.
        else:
            pathway_network = reactome_network.get_network(
                network.nodes(),
                enrichment_test=test.ENRICHMENT_TEST[(
                    configuration["Reactome network"].get(
                        "test", "hypergeometric"),
                    configuration["Reactome network"].get("increase", True))],
                multiple_testing_correction=correction.CORRECTION[
                    configuration["Reactome network"].get(
                        "correction", "Benjamini-Yekutieli")],
                organism=configuration["Reactome network"].get(
                    "organism", 9606),
                file_pathways=configuration.get("Reactome",
                                                {}).get("file",
                                                        {}).get("pathways"),
                file_pathways_relation=configuration.get("Reactome", {}).get(
                    "file", {}).get("pathways relation"),
                file_accession_map=configuration.get("Reactome", {}).get(
                    "file", {}).get("accession map"),
                file_uniprot=configuration.get("UniProt", {}).get("file"))

        # Export the Reactome network.
        file = reactome_network.export(pathway_network,
                                       f"{identifier}_reactome")

        if file is None:
            logger.warning(
                "The Reactome network was not exported due to naming conflict.")
        else:
            logger.info("The Reactome network was exported to %s.", file)

        # Export the Cytoscape style for the Reactome network.
        if "Cytoscape" in configuration:
            pathway_network_styles = reactome_network_style.get_styles(
                pathway_network)
            file = reactome_network_style.export(pathway_network_styles,
                                                 f"{identifier}_reactome")

            if file is None:
                logger.warning(
                    "The Cytoscape style for the Reactome network was not "
                    "exported due to naming conflict.")
            else:
                logger.info(
                    "The Cytoscape style for the Reactome network was exported "
                    "to %s.", file)

    # Analyze communities of the protein-protein interaction network.
    if "community detection" in configuration:
        logger.info("Detecting communities.")

        # Derive edge weights for community detection from protein-protein
        # interaction confidence scores of the protein-protein interaction
        # network.
        protein_interaction_network.set_edge_weights(
            network,
            weight=average.CONFIDENCE_SCORE_AVERAGE[
                configuration["community detection"].get("edge weight")])

        # Identify communities of the protein-protein interaction network.
        communities = protein_interaction_network.get_communities(
            network,
            community_size=configuration["community detection"].get(
                "community size", network.number_of_nodes()),
            community_size_average=average.MODULE_SIZE_AVERAGE[
                configuration["community detection"].get(
                    "community size average", "mean")],
            algorithm=modularization.ALGORITHM[
                configuration["community detection"].get(
                    "algorithm", "Louvain")],
            resolution=configuration["community detection"].get(
                "resolution", 1.0))

        protein_interaction_network.remove_edge_weights(network)

        export = {community: False for community in communities}

    # Assess Gene Ontology term enrichment by the protein-protein interaction
    # network or its densely interacting communities.
    if (("Gene Ontology enrichment" in configuration or
         "Gene Ontology enrichment" in configuration.get(
             "community detection", {})) and
            not os.path.isfile(f"{identifier}_gene_ontology.tsv")):
        logger.info("Gene Ontology enrichment test results are exported to %s.",
                    f"{identifier}_gene_ontology.tsv")

        with open(f"{identifier}_gene_ontology.tsv",
                  "w",
                  newline="",
                  encoding="utf-8") as gene_ontology_table:
            gene_ontology_writer = csv.writer(gene_ontology_table,
                                              dialect="excel-tab")
            gene_ontology_writer.writerow([
                "network", "term", "name", "p-value", "number of proteins",
                "number of associated proteins", "associated proteins"
            ])

            # Assess Gene Ontology term enrichment by subsets of proteins from
            # the protein-protein interaction network derived from average
            # measurements of the proteins for different types of
            # post-translational modification at different times of measurement.
            if "post-translational modifications" in configuration.get(
                    "Gene Ontology enrichment", {}):
                if configuration["Gene Ontology enrichment"].get(
                        "intersection", False):
                    proteins = set(network.nodes())
                else:
                    proteins = set()

                for time in protein_interaction_network.get_times(network):
                    for m in configuration["Gene Ontology enrichment"].get(
                            "post-translational modifications", []):
                        if m in protein_interaction_network.get_modifications(
                                network, time):
                            measurement_range = (
                                score.MEASUREMENT_SCORE[configuration[
                                    "Gene Ontology enrichment"].get(
                                        "score", {}).get(m)]
                                (configuration["Gene Ontology enrichment"].get(
                                    "measurement", {}).get(
                                        m,
                                        default.MEASUREMENT_RANGE[configuration[
                                            "Gene Ontology enrichment"].get(
                                                "score", {}).get(m)])[0],
                                 protein_interaction_network.get_measurements(
                                     network, time, m,
                                     average.SITE_AVERAGE[configuration[
                                         "Gene Ontology enrichment"].get(
                                             "site average", {}).get(
                                                 m,
                                                 "maximum absolute logarithm")],
                                     average.REPLICATE_AVERAGE[configuration[
                                         "Gene Ontology enrichment"].get(
                                             "replicate average",
                                             {}).get(m, "mean")])),
                                score.MEASUREMENT_SCORE[configuration[
                                    "Gene Ontology enrichment"].get(
                                        "score", {}).get(m)]
                                (configuration["Gene Ontology enrichment"].get(
                                    "measurement", {}).get(
                                        m,
                                        default.MEASUREMENT_RANGE[configuration[
                                            "Gene Ontology enrichment"].get(
                                                "score", {}).get(m)])[1],
                                 protein_interaction_network.get_measurements(
                                     network, time, m,
                                     average.SITE_AVERAGE[configuration[
                                         "Gene Ontology enrichment"].get(
                                             "site average", {}).get(
                                                 m,
                                                 "maximum absolute logarithm")],
                                     average.REPLICATE_AVERAGE[configuration[
                                         "Gene Ontology enrichment"].get(
                                             "replicate average",
                                             {}).get(m, "mean")])))

                            if configuration["Gene Ontology enrichment"].get(
                                    "intersection", False):
                                proteins.intersection_update(
                                    protein_interaction_network.get_proteins(
                                        network, time, m,
                                        average.SITE_AVERAGE[configuration[
                                            "Gene Ontology enrichment"].get(
                                                "site average", {}).get(
                                                    m,
                                                    "maximum absolute logarithm"
                                                )],
                                        average.REPLICATE_AVERAGE[configuration[
                                            "Gene Ontology enrichment"].get(
                                                "replicate average",
                                                {}).get(m, "mean")],
                                        measurement_range))

                            else:
                                proteins.update(
                                    protein_interaction_network.get_proteins(
                                        network, time, m,
                                        average.SITE_AVERAGE[configuration[
                                            "Gene Ontology enrichment"].get(
                                                "site average", {}).get(
                                                    m,
                                                    "maximum absolute logarithm"
                                                )],
                                        average.REPLICATE_AVERAGE[configuration[
                                            "Gene Ontology enrichment"].get(
                                                "replicate average",
                                                {}).get(m, "mean")],
                                        measurement_range))

                gene_ontology_enrichment = gene_ontology.get_enrichment(
                    [frozenset(proteins)],
                    reference=[set()]
                    if configuration["Gene Ontology enrichment"].get(
                        "annotation", False) else [network.nodes()],
                    enrichment_test=test.ENRICHMENT_TEST[(
                        configuration["Gene Ontology enrichment"].get(
                            "test", "hypergeometric"),
                        configuration["Gene Ontology enrichment"].get(
                            "increase", True))],
                    multiple_testing_correction=correction.CORRECTION[
                        configuration["Gene Ontology enrichment"].get(
                            "correction", "Benjamini-Yekutieli")],
                    organism=configuration["Gene Ontology enrichment"].get(
                        "organism", 9606),
                    namespaces=[
                        namespace.replace(" ", "_") for namespace in
                        configuration["Gene Ontology enrichment"].get(
                            "namespaces", [])
                    ],
                    file_ontology=configuration.get("Gene Ontology",
                                                    {}).get("file",
                                                            {}).get("ontology"),
                    file_annotation=configuration.get("Gene Ontology", {}).get(
                        "file", {}).get("annotation"),
                    file_annotation_isoform=configuration.get(
                        "Gene Ontology", {}).get("file",
                                                 {}).get("annotation isoform"),
                    file_uniprot=configuration.get("UniProt", {}).get("file"))

                for (term, name), (p, prt) in sorted(
                        gene_ontology_enrichment[frozenset(proteins)].items(),
                        key=lambda item: item[0][0]):
                    if p <= configuration["Gene Ontology enrichment"].get(
                            "p", 1.0):
                        gene_ontology_writer.writerow([
                            "network", term, name, p,
                            network.number_of_nodes(),
                            len(prt), " ".join(sorted(prt))
                        ])

            # Assess Gene Ontology term enrichment by the protein-protein
            # interaction network.
            elif "Gene Ontology enrichment" in configuration:
                gene_ontology_enrichment = gene_ontology.get_enrichment(
                    [frozenset(network.nodes())],
                    reference=[set()],
                    enrichment_test=test.ENRICHMENT_TEST[(
                        configuration["Gene Ontology enrichment"].get(
                            "test", "hypergeometric"),
                        configuration["Gene Ontology enrichment"].get(
                            "increase", True))],
                    multiple_testing_correction=correction.CORRECTION[
                        configuration["Gene Ontology enrichment"].get(
                            "correction", "Benjamini-Yekutieli")],
                    organism=configuration["Gene Ontology enrichment"].get(
                        "organism", 9606),
                    namespaces=[
                        namespace.replace(" ", "_") for namespace in
                        configuration["Gene Ontology enrichment"].get(
                            "namespaces", [])
                    ],
                    file_ontology=configuration.get("Gene Ontology",
                                                    {}).get("file",
                                                            {}).get("ontology"),
                    file_annotation=configuration.get("Gene Ontology", {}).get(
                        "file", {}).get("annotation"),
                    file_annotation_isoform=configuration.get(
                        "Gene Ontology", {}).get("file",
                                                 {}).get("annotation isoform"),
                    file_uniprot=configuration.get("UniProt", {}).get("file"))

                for (term,
                     name), (p,
                             prt) in sorted(gene_ontology_enrichment[frozenset(
                                 network.nodes())].items(),
                                            key=lambda item: item[0][0]):
                    if p <= configuration["Gene Ontology enrichment"].get(
                            "p", 1.0):
                        gene_ontology_writer.writerow([
                            "network", term, name, p,
                            network.number_of_nodes(),
                            len(prt), " ".join(sorted(prt))
                        ])

            # Assess local Gene Ontology term enrichment by subsets of proteins
            # from communities of the protein-protein interaction network
            # derived from average measurements of the proteins for different
            # types of post-translational modification at different times of
            # measurement.
            if "post-translational modifications" in configuration.get(
                    "community detection", {}).get("Gene Ontology enrichment",
                                                   {}):
                if configuration["community detection"][
                        "Gene Ontology enrichment"].get("intersection", False):
                    subsets = {
                        community: set(community.nodes())
                        for community in communities
                    }
                else:
                    subsets = {community: set() for community in communities}

                for time in protein_interaction_network.get_times(network):
                    for m in configuration["community detection"][
                            "Gene Ontology enrichment"].get(
                                "post-translational modifications", []):
                        if m in protein_interaction_network.get_modifications(
                                network, time):
                            for community in subsets:
                                measurement_range = (
                                    score.MEASUREMENT_SCORE[
                                        configuration["community detection"]
                                        ["Gene Ontology enrichment"].get(
                                            "score", {}).get(m)]
                                    (configuration["community detection"]
                                     ["Gene Ontology enrichment"].get(
                                         "measurement", {}).get(
                                             m, default.MEASUREMENT_RANGE[
                                                 configuration[
                                                     "community detection"]
                                                 ["Gene Ontology enrichment"].
                                                 get("score", {}).get(m)])[0],
                                     protein_interaction_network.
                                     get_measurements(
                                         network, time, m, average.SITE_AVERAGE[
                                             configuration[
                                                 "community detection"]
                                             ["Gene Ontology enrichment"].get(
                                                 "site average", {}).
                                             get(m,
                                                 "maximum absolute logarithm")],
                                         average.REPLICATE_AVERAGE[
                                             configuration[
                                                 "community detection"]
                                             ["Gene Ontology enrichment"].get(
                                                 "replicate average",
                                                 {}).get(m, "mean")])),
                                    score.MEASUREMENT_SCORE[
                                        configuration["community detection"]
                                        ["Gene Ontology enrichment"].get(
                                            "score", {}).get(m)]
                                    (configuration["community detection"]
                                     ["Gene Ontology enrichment"].get(
                                         "measurement", {}).get(
                                             m, default.MEASUREMENT_RANGE[
                                                 configuration[
                                                     "community detection"]
                                                 ["Gene Ontology enrichment"].
                                                 get("score", {}).get(m)])[1],
                                     protein_interaction_network.
                                     get_measurements(
                                         network, time, m, average.SITE_AVERAGE[
                                             configuration[
                                                 "community detection"]
                                             ["Gene Ontology enrichment"].get(
                                                 "site average", {}).
                                             get(m,
                                                 "maximum absolute logarithm")],
                                         average.REPLICATE_AVERAGE[
                                             configuration[
                                                 "community detection"]
                                             ["Gene Ontology enrichment"].get(
                                                 "replicate average",
                                                 {}).get(m, "mean")])))

                                if configuration["community detection"][
                                        "Gene Ontology enrichment"].get(
                                            "intersection", False):
                                    subsets[community].intersection_update(
                                        protein_interaction_network.
                                        get_proteins(
                                            community, time, m,
                                            average.SITE_AVERAGE[
                                                configuration[
                                                    "community detection"]
                                                ["Gene Ontology enrichment"].
                                                get("site average", {}).get(
                                                    m,
                                                    "maximum absolute logarithm"
                                                )],
                                            average.REPLICATE_AVERAGE[
                                                configuration[
                                                    "community detection"]
                                                ["Gene Ontology enrichment"].
                                                get("replicate average",
                                                    {}).get(m, "mean")],
                                            measurement_range))

                                else:
                                    subsets[community].update(
                                        protein_interaction_network.
                                        get_proteins(
                                            community, time, m,
                                            average.SITE_AVERAGE[
                                                configuration[
                                                    "community detection"]
                                                ["Gene Ontology enrichment"].
                                                get("site average", {}).get(
                                                    m,
                                                    "maximum absolute logarithm"
                                                )],
                                            average.REPLICATE_AVERAGE[
                                                configuration[
                                                    "community detection"]
                                                ["Gene Ontology enrichment"].
                                                get("replicate average",
                                                    {}).get(m, "mean")],
                                            measurement_range))

                gene_ontology_enrichment = gene_ontology.get_enrichment(
                    [
                        frozenset(subsets[community])
                        for community in communities
                    ],
                    reference=[network.nodes()]
                    if configuration["community detection"]
                    ["Gene Ontology enrichment"].get("network", False) else
                    ([set()] if configuration["community detection"]
                     ["Gene Ontology enrichment"].get("annotation", False) else
                     [community.nodes() for community in communities]),
                    enrichment_test=test.ENRICHMENT_TEST[(
                        configuration["community detection"]
                        ["Gene Ontology enrichment"].get(
                            "test", "hypergeometric"),
                        configuration["community detection"]
                        ["Gene Ontology enrichment"].get("increase", True))],
                    multiple_testing_correction=correction.CORRECTION[
                        configuration["community detection"]
                        ["Gene Ontology enrichment"].get(
                            "correction", "Benjamini-Yekutieli")],
                    organism=configuration["community detection"]
                    ["Gene Ontology enrichment"].get("organism", 9606),
                    namespaces=[
                        namespace.replace(" ", "_")
                        for namespace in configuration["community detection"]
                        ["Gene Ontology enrichment"].get("namespaces", [])
                    ],
                    file_ontology=configuration.get("Gene Ontology",
                                                    {}).get("file",
                                                            {}).get("ontology"),
                    file_annotation=configuration.get("Gene Ontology", {}).get(
                        "file", {}).get("annotation"),
                    file_annotation_isoform=configuration.get(
                        "Gene Ontology", {}).get("file",
                                                 {}).get("annotation isoform"),
                    file_uniprot=configuration.get("UniProt", {}).get("file"))

                for k, community in enumerate(sorted(
                        communities,
                        key=lambda community: int(community.number_of_nodes()),
                        reverse=True),
                                              start=1):
                    for (term, name), (p, prt) in sorted(
                            gene_ontology_enrichment[frozenset(
                                subsets[community])].items(),
                            key=lambda item: item[0][0]):
                        if p <= configuration["community detection"][
                                "Gene Ontology enrichment"].get("p", 1.0):
                            export[community] = True
                            gene_ontology_writer.writerow([
                                f"community {k}", term, name, p,
                                community.number_of_nodes(),
                                len(prt), " ".join(sorted(prt))
                            ])

            # Assess local Gene Ontology term enrichment by communities of the
            # protein-protein interaction network.
            elif "Gene Ontology enrichment" in configuration.get(
                    "community detection", {}):
                gene_ontology_enrichment = gene_ontology.get_enrichment(
                    [frozenset(community.nodes()) for community in communities],
                    reference=[set()] if configuration["community detection"]
                    ["Gene Ontology enrichment"].get(
                        "annotation", False) else [network.nodes()],
                    enrichment_test=test.ENRICHMENT_TEST[(
                        configuration["community detection"]
                        ["Gene Ontology enrichment"].get(
                            "test", "hypergeometric"),
                        configuration["community detection"]
                        ["Gene Ontology enrichment"].get("increase", True))],
                    multiple_testing_correction=correction.CORRECTION[
                        configuration["community detection"]
                        ["Gene Ontology enrichment"].get(
                            "correction", "Benjamini-Yekutieli")],
                    organism=configuration["community detection"]
                    ["Gene Ontology enrichment"].get("organism", 9606),
                    namespaces=[
                        namespace.replace(" ", "_")
                        for namespace in configuration["community detection"]
                        ["Gene Ontology enrichment"].get("namespaces", [])
                    ],
                    file_ontology=configuration.get("Gene Ontology",
                                                    {}).get("file",
                                                            {}).get("ontology"),
                    file_annotation=configuration.get("Gene Ontology", {}).get(
                        "file", {}).get("annotation"),
                    file_annotation_isoform=configuration.get(
                        "Gene Ontology", {}).get("file",
                                                 {}).get("annotation isoform"),
                    file_uniprot=configuration.get("UniProt", {}).get("file"))

                for k, community in enumerate(sorted(
                        communities,
                        key=lambda community: int(community.number_of_nodes()),
                        reverse=True),
                                              start=1):
                    for (term, name), (p, prt) in sorted(
                            gene_ontology_enrichment[frozenset(
                                community.nodes())].items(),
                            key=lambda item: item[0][0]):
                        if p <= configuration["community detection"][
                                "Gene Ontology enrichment"].get("p", 1.0):
                            export[community] = True
                            gene_ontology_writer.writerow([
                                f"community {k}", term, name, p,
                                community.number_of_nodes(),
                                len(prt), " ".join(sorted(prt))
                            ])
    else:
        logger.warning(
            "Gene Ontology enrichment test results are not exported due to "
            "naming conflict.")

    # Assess Reactome pathway enrichment by the protein-protein interaction
    # network or its densely interacting communities.
    if (("Reactome enrichment" in configuration or "Reactome enrichment"
         in configuration.get("community detection", {})) and
            not os.path.isfile(f"{identifier}_reactome.tsv")):
        logger.info("Reactome enrichment test results are exported to %s.",
                    f"{identifier}_reactome.tsv")

        with open(f"{identifier}_reactome.tsv",
                  "w",
                  newline="",
                  encoding="utf-8") as reactome_table:
            reactome_writer = csv.writer(reactome_table, dialect="excel-tab")
            reactome_writer.writerow([
                "network", "pathway", "name", "p-value", "number of proteins",
                "number of associated proteins", "associated proteins"
            ])

            # Assess Reactome pathway enrichment by subsets of proteins from the
            # protein-protein interaction network derived from average
            # measurements of the proteins for different types of
            # post-translational modification at different times of measurement.
            if "post-translational modifications" in configuration.get(
                    "Reactome enrichment", {}):
                if configuration["Reactome enrichment"].get(
                        "intersection", False):
                    proteins = set(network.nodes())
                else:
                    proteins = set()

                for time in protein_interaction_network.get_times(network):
                    for m in configuration["Reactome enrichment"].get(
                            "post-translational modifications", []):
                        if m in protein_interaction_network.get_modifications(
                                network, time):
                            measurement_range = (
                                score.MEASUREMENT_SCORE[
                                    configuration["Reactome enrichment"].get(
                                        "score", {}).get(m)]
                                (configuration["Reactome enrichment"].get(
                                    "measurement", {}).get(
                                        m,
                                        default.MEASUREMENT_RANGE[configuration[
                                            "Reactome enrichment"].get(
                                                "score", {}).get(m)])[0],
                                 protein_interaction_network.get_measurements(
                                     network, time, m,
                                     average.SITE_AVERAGE[configuration[
                                         "Reactome enrichment"].get(
                                             "site average", {}).get(
                                                 m,
                                                 "maximum absolute logarithm")],
                                     average.REPLICATE_AVERAGE[configuration[
                                         "Reactome enrichment"].get(
                                             "replicate average",
                                             {}).get(m, "mean")])),
                                score.MEASUREMENT_SCORE[
                                    configuration["Reactome enrichment"].get(
                                        "score", {}).get(m)]
                                (configuration["Reactome enrichment"].get(
                                    "measurement", {}).get(
                                        m,
                                        default.MEASUREMENT_RANGE[configuration[
                                            "Reactome enrichment"].get(
                                                "score", {}).get(m)])[1],
                                 protein_interaction_network.get_measurements(
                                     network, time, m,
                                     average.SITE_AVERAGE[configuration[
                                         "Reactome enrichment"].get(
                                             "site average", {}).get(
                                                 m,
                                                 "maximum absolute logarithm")],
                                     average.REPLICATE_AVERAGE[configuration[
                                         "Reactome enrichment"].get(
                                             "replicate average",
                                             {}).get(m, "mean")])))

                            if configuration["Reactome enrichment"].get(
                                    "intersection", False):
                                proteins.intersection_update(
                                    protein_interaction_network.get_proteins(
                                        network, time, m, average.SITE_AVERAGE[
                                            configuration["Reactome enrichment"]
                                            .get("site average", {}).get(
                                                m,
                                                "maximum absolute logarithm")],
                                        average.REPLICATE_AVERAGE[configuration[
                                            "Reactome enrichment"].get(
                                                "replicate average",
                                                {}).get(m, "mean")],
                                        measurement_range))

                            else:
                                proteins.update(
                                    protein_interaction_network.get_proteins(
                                        network, time, m, average.SITE_AVERAGE[
                                            configuration["Reactome enrichment"]
                                            .get("site average", {}).get(
                                                m,
                                                "maximum absolute logarithm")],
                                        average.REPLICATE_AVERAGE[configuration[
                                            "Reactome enrichment"].get(
                                                "replicate average",
                                                {}).get(m, "mean")],
                                        measurement_range))

                reactome_enrichment = reactome.get_enrichment(
                    [frozenset(proteins)],
                    reference=[set()]
                    if configuration["Reactome enrichment"].get(
                        "annotation", False) else [network.nodes()],
                    enrichment_test=test.ENRICHMENT_TEST[(
                        configuration["Reactome enrichment"].get(
                            "test", "hypergeometric"),
                        configuration["Reactome enrichment"].get(
                            "increase", True))],
                    multiple_testing_correction=correction.CORRECTION[
                        configuration["Reactome enrichment"].get(
                            "correction", "Benjamini-Yekutieli")],
                    organism=configuration["Reactome enrichment"].get(
                        "organism", 9606),
                    file_pathways=configuration.get("Reactome",
                                                    {}).get("file",
                                                            {}).get("pathways"),
                    file_accession_map=configuration.get("Reactome", {}).get(
                        "file", {}).get("accession map"),
                    file_uniprot=configuration.get("UniProt", {}).get("file"))

                for (pathway, name), (p, prt) in sorted(
                        reactome_enrichment[frozenset(proteins)].items(),
                        key=lambda item: item[0][0]):
                    if p <= configuration["Reactome enrichment"].get("p", 1.0):
                        reactome_writer.writerow([
                            "network", pathway, name, p,
                            network.number_of_nodes(),
                            len(prt), " ".join(sorted(prt))
                        ])

            # Assess Reactome pathway enrichment by the protein-protein
            # interaction network.
            elif "Reactome enrichment" in configuration:
                reactome_enrichment = reactome.get_enrichment(
                    [frozenset(network.nodes())],
                    reference=[set()],
                    enrichment_test=test.ENRICHMENT_TEST[(
                        configuration["Reactome enrichment"].get(
                            "test", "hypergeometric"),
                        configuration["Reactome enrichment"].get(
                            "increase", True))],
                    multiple_testing_correction=correction.CORRECTION[
                        configuration["Reactome enrichment"].get(
                            "correction", "Benjamini-Yekutieli")],
                    organism=configuration["Reactome enrichment"].get(
                        "organism", 9606),
                    file_pathways=configuration.get("Reactome",
                                                    {}).get("file",
                                                            {}).get("pathways"),
                    file_accession_map=configuration.get("Reactome", {}).get(
                        "file", {}).get("accession map"),
                    file_uniprot=configuration.get("UniProt", {}).get("file"))

                for (pathway, name), (p, prt) in sorted(
                        reactome_enrichment[frozenset(network.nodes())].items(),
                        key=lambda item: item[0][0]):
                    if p <= configuration["Reactome enrichment"].get("p", 1.0):
                        reactome_writer.writerow([
                            "network", pathway, name, p,
                            network.number_of_nodes(),
                            len(prt), " ".join(sorted(prt))
                        ])

            # Assess local Reactome pathway enrichment by subsets of proteins
            # from communities of the protein-protein interaction network
            # derived from average measurements of the proteins for different
            # types of post-translational modification at different times of
            # measurement.
            if "post-translational modifications" in configuration.get(
                    "community detection", {}).get("Reactome enrichment", {}):
                if configuration["community detection"][
                        "Reactome enrichment"].get("intersection", False):
                    subsets = {
                        community: set(community.nodes())
                        for community in communities
                    }
                else:
                    subsets = {community: set() for community in communities}

                for time in protein_interaction_network.get_times(network):
                    for m in configuration["community detection"][
                            "Reactome enrichment"].get(
                                "post-translational modifications", []):
                        if m in protein_interaction_network.get_modifications(
                                network, time):
                            for community in subsets:
                                measurement_range = (
                                    score.MEASUREMENT_SCORE[
                                        configuration["community detection"]
                                        ["Reactome enrichment"].get(
                                            "score", {}).get(m)]
                                    (configuration["community detection"]
                                     ["Reactome enrichment"].get(
                                         "measurement", {}).get(
                                             m, default.MEASUREMENT_RANGE[
                                                 configuration[
                                                     "community detection"]
                                                 ["Reactome enrichment"].get(
                                                     "score", {}).get(m)])[0],
                                     protein_interaction_network.
                                     get_measurements(
                                         network, time, m, average.SITE_AVERAGE[
                                             configuration[
                                                 "community detection"]
                                             ["Reactome enrichment"].get(
                                                 "site average", {}).
                                             get(m,
                                                 "maximum absolute logarithm")],
                                         average.REPLICATE_AVERAGE[
                                             configuration[
                                                 "community detection"]
                                             ["Reactome enrichment"].get(
                                                 "replicate average",
                                                 {}).get(m, "mean")])),
                                    score.MEASUREMENT_SCORE[
                                        configuration["community detection"]
                                        ["Reactome enrichment"].get(
                                            "score", {}).get(m)]
                                    (configuration["community detection"]
                                     ["Reactome enrichment"].get(
                                         "measurement", {}).get(
                                             m, default.MEASUREMENT_RANGE[
                                                 configuration[
                                                     "community detection"]
                                                 ["Reactome enrichment"].get(
                                                     "score", {}).get(m)])[1],
                                     protein_interaction_network.
                                     get_measurements(
                                         network, time, m, average.SITE_AVERAGE[
                                             configuration[
                                                 "community detection"]
                                             ["Reactome enrichment"].get(
                                                 "site average", {}).
                                             get(m,
                                                 "maximum absolute logarithm")],
                                         average.REPLICATE_AVERAGE[
                                             configuration[
                                                 "community detection"]
                                             ["Reactome enrichment"].get(
                                                 "replicate average",
                                                 {}).get(m, "mean")])))

                                if configuration["community detection"][
                                        "Reactome enrichment"].get(
                                            "intersection", False):
                                    subsets[community].intersection_update(
                                        protein_interaction_network.
                                        get_proteins(
                                            community, time, m,
                                            average.SITE_AVERAGE[
                                                configuration[
                                                    "community detection"]
                                                ["Reactome enrichment"].
                                                get("site average", {}).get(
                                                    m,
                                                    "maximum absolute logarithm"
                                                )], average.REPLICATE_AVERAGE[
                                                    configuration[
                                                        "community detection"]
                                                    ["Reactome enrichment"].get(
                                                        "replicate average",
                                                        {}).get(m, "mean")],
                                            measurement_range))

                                else:
                                    subsets[community].update(
                                        protein_interaction_network.
                                        get_proteins(
                                            community, time, m,
                                            average.SITE_AVERAGE[
                                                configuration[
                                                    "community detection"]
                                                ["Reactome enrichment"].
                                                get("site average", {}).get(
                                                    m,
                                                    "maximum absolute logarithm"
                                                )], average.REPLICATE_AVERAGE[
                                                    configuration[
                                                        "community detection"]
                                                    ["Reactome enrichment"].get(
                                                        "replicate average",
                                                        {}).get(m, "mean")],
                                            measurement_range))

                reactome_enrichment = reactome.get_enrichment(
                    [
                        frozenset(subsets[community])
                        for community in communities
                    ],
                    reference=[network.nodes()]
                    if configuration["community detection"]
                    ["Reactome enrichment"].get("network", False) else
                    ([set()] if configuration["community detection"]
                     ["Reactome enrichment"].get("annotation", False) else
                     [community.nodes() for community in communities]),
                    enrichment_test=test.ENRICHMENT_TEST[(
                        configuration["community detection"]
                        ["Reactome enrichment"].get("test", "hypergeometric"),
                        configuration["community detection"]
                        ["Reactome enrichment"].get("increase", True))],
                    multiple_testing_correction=correction.CORRECTION[
                        configuration["community detection"]
                        ["Reactome enrichment"].get("correction",
                                                    "Benjamini-Yekutieli")],
                    organism=configuration["community detection"]
                    ["Reactome enrichment"].get("organism", 9606),
                    file_pathways=configuration.get("Reactome",
                                                    {}).get("file",
                                                            {}).get("pathways"),
                    file_accession_map=configuration.get("Reactome", {}).get(
                        "file", {}).get("accession map"),
                    file_uniprot=configuration.get("UniProt", {}).get("file"))

                for k, community in enumerate(sorted(
                        communities,
                        key=lambda community: int(community.number_of_nodes()),
                        reverse=True),
                                              start=1):
                    for (pathway,
                         name), (p,
                                 prt) in sorted(reactome_enrichment[frozenset(
                                     subsets[community])].items(),
                                                key=lambda item: item[0][0]):
                        if p <= configuration["community detection"][
                                "Reactome enrichment"].get("p", 1.0):
                            export[community] = True
                            reactome_writer.writerow([
                                f"community {k}", pathway, name, p,
                                community.number_of_nodes(),
                                len(prt), " ".join(sorted(prt))
                            ])

            # Assess local Reactome pathway enrichment by communities of the
            # protein-protein interaction network.
            elif "Reactome enrichment" in configuration:
                reactome_enrichment = reactome.get_enrichment(
                    [frozenset(community.nodes()) for community in communities],
                    reference=[set()] if configuration["community detection"]
                    ["Reactome enrichment"].get("annotation",
                                                False) else [network.nodes()],
                    enrichment_test=test.ENRICHMENT_TEST[(
                        configuration["community detection"]
                        ["Reactome enrichment"].get("test", "hypergeometric"),
                        configuration["community detection"]
                        ["Reactome enrichment"].get("increase", True))],
                    multiple_testing_correction=correction.CORRECTION[
                        configuration["community detection"]
                        ["Reactome enrichment"].get("correction",
                                                    "Benjamini-Yekutieli")],
                    organism=configuration["community detection"]
                    ["Reactome enrichment"].get("organism", 9606),
                    file_pathways=configuration.get("Reactome",
                                                    {}).get("file",
                                                            {}).get("pathways"),
                    file_accession_map=configuration.get("Reactome", {}).get(
                        "file", {}).get("accession map"),
                    file_uniprot=configuration.get("UniProt", {}).get("file"))

                for k, community in enumerate(sorted(
                        communities,
                        key=lambda community: int(community.number_of_nodes()),
                        reverse=True),
                                              start=1):
                    for (pathway,
                         name), (p,
                                 prt) in sorted(reactome_enrichment[frozenset(
                                     community.nodes())].items(),
                                                key=lambda item: item[0][0]):
                        if p <= configuration["community detection"][
                                "Reactome enrichment"].get("p", 1.0):
                            export[community] = True
                            reactome_writer.writerow([
                                f"community {k}", pathway, name, p,
                                community.number_of_nodes(),
                                len(prt), " ".join(sorted(prt))
                            ])

    else:
        logger.warning(
            "Reactome enrichment test results are not exported due to naming "
            "conflict.")

    # Assess the enrichment of average measurements a at most a lower or at
    # least an upper threshold by communities of the protein-protein interaction
    # network for different types of post-translational modification at
    # different times of measurement.
    if ("measurement enrichment" in configuration.get("community detection", {})
            and not os.path.isfile(f"{identifier}_measurement_enrichment.tsv")):
        logger.info("Measurement enrichment test results are exported to %s.",
                    f"{identifier}_measurement_enrichment.tsv")

        with open(f"{identifier}_measurement_enrichment.tsv",
                  "w",
                  newline="",
                  encoding="utf-8") as measurement_enrichment_table:
            measurement_enrichment_writer = csv.writer(
                measurement_enrichment_table, dialect="excel-tab")
            measurement_enrichment_writer.writerow([
                "community", "time", "post-translational modification",
                "p-value", "number of proteins",
                "number of associated proteins", "associated proteins"
            ])

            measurement_ranges = {
                modification:
                default.MEASUREMENT_RANGE.get(score,
                                              default.MEASUREMENT_RANGE[None])
                for modification, score in configuration["community detection"]
                ["measurement enrichment"].get("score", {}).items()
            }

            for modification, measurement_range in configuration[
                    "community detection"]["measurement enrichment"].get(
                        "measurement", {}).items():
                measurement_ranges[modification] = measurement_range

            enrichment = protein_interaction_network.get_enrichment(
                network,
                communities,
                measurement_ranges=measurement_ranges,
                measurement_score={
                    modification:
                    score.MEASUREMENT_SCORE.get(measurement_score,
                                                score.MEASUREMENT_SCORE[None])
                    for modification, measurement_score in
                    configuration["community detection"]
                    ["measurement enrichment"].get("score", {})
                },
                site_average={
                    modification: average.SITE_AVERAGE.get(
                        site_average,
                        average.SITE_AVERAGE["maximum absolute logarithm"])
                    if site_average is not None else None for modification,
                    site_average in configuration["community detection"]
                    ["measurement enrichment"].get("site average", {}).items()
                },
                replicate_average={
                    modification: average.REPLICATE_AVERAGE.get(
                        replicate_average, average.REPLICATE_AVERAGE["mean"])
                    if replicate_average is not None else None for modification,
                    replicate_average in configuration["community detection"]
                    ["measurement enrichment"].get("replicate average",
                                                   {}).items()
                },
                enrichment_test=test.ENRICHMENT_TEST[(
                    configuration["community detection"]
                    ["measurement enrichment"].get("test", "hypergeometric"),
                    configuration["community detection"]
                    ["measurement enrichment"].get("increase", True))],
                multiple_testing_correction=correction.CORRECTION[
                    configuration["community detection"]
                    ["measurement enrichment"].get("correction",
                                                   "Benjamini-Yekutieli")])

            for k, community in enumerate(sorted(
                    communities,
                    key=lambda community: int(community.number_of_nodes()),
                    reverse=True),
                                          start=1):
                for time in enrichment[community]:
                    for modification, (p, prt) in sorted(
                            enrichment[community][time].items(),
                            key=lambda item: item[0][0]):
                        if p <= configuration["community detection"][
                                "measurement enrichment"].get("p", 1.0):
                            export[community] = True
                            measurement_enrichment_writer.writerow([
                                f"community {k}", time, modification, p,
                                community.number_of_nodes(),
                                len(prt), " ".join(sorted(prt))
                            ])
    else:
        logger.warning(
            "Measurement enrichment test results are not exported due to "
            "naming conflict.")

    # Assess the equality of location of average measurements of communities
    # relative to the remaining protein-protein interaction network for
    # different types of post-translational modification at different times of
    # measurement.
    if ("measurement location" in configuration.get("community detection", {})
            and not os.path.isfile(f"{identifier}_measurement_location.tsv")):
        logger.info("Measurement location test results are exported to %s.",
                    f"{identifier}_measurement_location.tsv")

        with open(f"{identifier}_measurement_location.tsv",
                  "w",
                  newline="",
                  encoding="utf-8") as measurement_location_table:
            measurement_location_writer = csv.writer(measurement_location_table,
                                                     dialect="excel-tab")
            measurement_location_writer.writerow([
                "community", "time", "post-translational modification",
                "p-value", "number of proteins"
            ])

            location = protein_interaction_network.get_location(
                network,
                communities,
                site_average={
                    modification: average.SITE_AVERAGE.get(
                        site_average,
                        average.SITE_AVERAGE["maximum absolute logarithm"])
                    if site_average is not None else None for modification,
                    site_average in configuration["community detection"]
                    ["measurement location"].get("site average", {}).items()
                },
                replicate_average={
                    modification: average.REPLICATE_AVERAGE.get(
                        replicate_average, average.REPLICATE_AVERAGE["mean"])
                    if replicate_average is not None else None
                    for modification, replicate_average in
                    configuration["community detection"]["measurement location"]
                    .get("replicate average", {}).items()
                },
                location_test=test.LOCATION_TEST[(
                    configuration["community detection"]
                    ["measurement location"].get("test",
                                                 "Mann-Whitney-Wilcoxon"),
                    configuration["community detection"]
                    ["measurement location"].get("increase", True),
                    configuration["community detection"]
                    ["measurement location"].get("absolute", True))],
                multiple_testing_correction=correction.CORRECTION[
                    configuration["community detection"]
                    ["measurement location"].get("correction",
                                                 "Benjamini-Yekutieli")])

            for k, community in enumerate(sorted(
                    communities,
                    key=lambda community: int(community.number_of_nodes()),
                    reverse=True),
                                          start=1):
                for time in location[community]:
                    for modification, p in sorted(
                            location[community][time].items(),
                            key=lambda item: item[0][0]):
                        if p <= configuration["community detection"][
                                "measurement location"].get("p", 1.0):
                            export[community] = True
                            measurement_location_writer.writerow([
                                f"community {k}", time, modification, p,
                                community.number_of_nodes()
                            ])

    else:
        logger.warning(
            "Measurement location test results are not exported due to naming "
            "conflict.")

    # Export the communities of the protein-protein interaction network
    # significant according to any of the specified hypothesis tests.
    files = []
    for k, community in enumerate(sorted(
            communities,
            key=lambda community: int(community.number_of_nodes()),
            reverse=True),
                                  start=1):
        if export[community]:
            files.append(
                protein_interaction_network.export(community,
                                                   f"{identifier}_{k}"))

    if any(export):
        if any(file is None for file in files):
            logger.warning(
                "Communities of the protein-protein interaction network were "
                "not exported due to naming conflicts.")
        else:
            logger.info(
                "Communities of the protein-protein interaction network were "
                "exported.")


def process_configuration(
        configurations: Mapping[str, Mapping[str, Any]]) -> None:
    """
    Executes workflows specified in configurations sequentially.

    Args:
        configurations: The specification of workflows.
    """
    # Process specified workflows sequentially.
    for identifier, configuration in configurations.items():
        process_workflow(identifier, configuration)


def process_configuration_file(configuration_file: str) -> None:
    """
    Launches execution of a workflow.

    Args:
        configuration_file: file name of configuration file.
    """
    # Process configuration file.
    with open(configuration_file, encoding="utf-8") as configuration:
        process_configuration(json.load(configuration))


def main() -> None:
    """Concurrent workflow execution."""
    # Specify parser for command line arguments.
    parser = argparse.ArgumentParser(
        description="DIANA: data integration and network analysis of "
        "post-translational modification based on mass spectrometry data")

    parser.add_argument("-c",
                        "--configuration",
                        help="configuration files",
                        nargs="+",
                        required=True)

    parser.add_argument(
        "-p",
        "--processes",
        help="maximum number of concurrently processed configuration files "
        f"(default: {os.cpu_count()})",
        type=int,
        default=os.cpu_count())

    parser.add_argument(
        "-l",
        "--logging",
        help="logging level (default: INFO)",
        default=logging.INFO,
        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"])

    args = parser.parse_args()

    # Configure logging.
    logging.basicConfig(
        filename=f"{os.path.splitext(os.path.basename(sys.argv[0]))[0]}.log",
        filemode="w",
        format="%(asctime)s %(levelname)s %(name)s (PID: %(process)d): "
        "%(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        level=args.logging,
        encoding="utf-8")

    # Process configuration files concurrently.
    with concurrent.futures.ProcessPoolExecutor(
            max_workers=args.processes) as executor:
        for process in executor.map(process_configuration_file,
                                    args.configuration):
            if process:
                process.results()


if __name__ == "__main__":
    main()
