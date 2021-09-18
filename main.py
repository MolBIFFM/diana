import argparse
import json
import logging
import os
import tempfile
import uuid

from pipeline.cytoscape_styles import CytoscapeStyles
from pipeline.databases import biogrid, corum, intact, reactome, string
from pipeline.interface import algorithm, combine, convert, extract
from pipeline.protein_interaction_network import ProteinInteractionNetwork
from pipeline.utilities import download, export


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-c",
        "--configurations",
        help="JSON configuration files",
        nargs="+",
        required=True,
    )
    parser.add_argument("--chunk-size",
                        help="file processing rate (default: {})".format(
                            download.CHUNK_SIZE),
                        type=int)
    parser.add_argument("-l", "--log", help="log file", type=str)
    parser.add_argument(
        "-u",
        "--user-agent",
        help="identity to download sources (default: {})".format(
            download.USER_AGENT),
        type=str)
    args = parser.parse_args()

    if args.log:
        logger = logging.getLogger("root")
        logging.basicConfig(filename=args.log,
                            filemode="w",
                            level=logging.DEBUG,
                            format="%(message)s")

    if args.chunk_size:
        download.CHUNK_SIZE = args.chunk_size

    if args.user_agent:
        download.USER_AGENT = args.user_agent

    for configuration_file in args.configurations:
        with open(configuration_file) as configuration:
            configurations = json.load(configuration)

        for i, configuration in enumerate(configurations, start=1):
            network = ProteinInteractionNetwork()

            for entry in configuration.get("genes", {}):
                if "file" in entry and "accession column" in entry:
                    network.add_genes_from_table(
                        file_name=entry["file"],
                        gene_accession_column=entry["accession column"],
                        gene_accession_format=extract.EXTRACT.get(
                            entry.get("accession format"),
                            lambda entry: [entry]),
                        sheet_name=entry.get("sheet", 0),
                        header=entry.get("header", 1) - 1,
                        taxon_identifier=entry.get("taxon identifier", 9606),
                    )

                elif "accessions" in entry:
                    network.add_genes_from(
                        genes=entry["accessions"],
                        taxon_identifier=entry.get("taxon identifier", 9606),
                    )

            for entry in configuration.get("proteins", {}):
                if "file" in entry and "accession column" in entry:
                    network.add_proteins_from_table(
                        file_name=entry["file"],
                        protein_accession_column=entry["accession column"],
                        protein_accession_format=extract.EXTRACT.get(
                            entry.get("accession format"),
                            lambda entry: [entry]),
                        time=entry.get("time", 0),
                        ptm=entry.get("post-translational modification", ""),
                        position_column=entry.get("position column", ""),
                        position_format=extract.EXTRACT.get(
                            entry.get("position format"),
                            lambda entry: [entry]),
                        replicates=entry.get("replicate columns", []),
                        sheet_name=entry.get("sheet", 0),
                        header=entry.get("header", 1) - 1,
                        num_sites=entry.get("sites", 1000),
                        num_replicates=entry.get("replicates", 1),
                        combine_replicates=combine.COMBINE_REPLICATES.get(
                            entry.get("combine replicates", "mean"),
                            combine.COMBINE_REPLICATES["mean"],
                        ),
                        convert_measurement=convert.LOG_BASE[entry.get(
                            "log base")],
                    )

                elif "accessions" in entry:
                    network.add_proteins_from(proteins=entry["accessions"], )

            if "protein-protein interactions" in configuration:
                k = 0
                while any(
                        configuration["protein-protein interactions"].get(
                            database, {}).get("neighbors", 0) > k
                        for database in
                    {"BioGRID", "CORUM", "IntAct", "Reactome", "STRING"}):
                    if "BioGRID" in configuration[
                            "protein-protein interactions"] and configuration[
                                "protein-protein interactions"]["BioGRID"].get(
                                    "neighbors", 0) > k:
                        biogrid.add_proteins(
                            network,
                            experimental_system=configuration[
                                "protein-protein interactions"]["BioGRID"].get(
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
                            experimental_system_type=configuration[
                                "protein-protein interactions"]["BioGRID"].get(
                                    "experimental system type",
                                    [
                                        "physical",
                                    ],
                                ),
                            taxon_identifier=configuration[
                                "protein-protein interactions"]["BioGRID"].get(
                                    "taxon identifier", 9606),
                            multi_validated_physical=configuration[
                                "protein-protein interactions"]["BioGRID"].get(
                                    "multi-validated physical", False),
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
                            mi_score=configuration[
                                "protein-protein interactions"]["IntAct"].get(
                                    "MI score", 0.27),
                        )

                    if "Reactome" in configuration[
                            "protein-protein interactions"] and configuration[
                                "protein-protein interactions"][
                                    "Reactome"].get("neighbors", 0) > k:
                        reactome.add_proteins(
                            network,
                            interaction_context=configuration[
                                "protein-protein interactions"]
                            ["Reactome"].get("interaction context", []),
                            interaction_type=configuration[
                                "protein-protein interactions"]
                            ["Reactome"].get("interaction type", []),
                            taxon_identifier=configuration[
                                "protein-protein interactions"]
                            ["Reactome"].get("taxon identifier", 9606),
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
                            fusion=configuration[
                                "protein-protein interactions"]["STRING"].get(
                                    "fusion", 0.0),
                            cooccurence=configuration[
                                "protein-protein interactions"]["STRING"].get(
                                    "cooccurence", 0.0),
                            homology=configuration[
                                "protein-protein interactions"]["STRING"].get(
                                    "homology", 0.0),
                            coexpression=configuration[
                                "protein-protein interactions"]["STRING"].get(
                                    "coexpression", 0.0),
                            coexpression_transferred=configuration[
                                "protein-protein interactions"]["STRING"].get(
                                    "coexpression transferred", 0.0),
                            experiments=configuration[
                                "protein-protein interactions"]["STRING"].get(
                                    "experiments", 0.7),
                            experiments_transferred=configuration[
                                "protein-protein interactions"]["STRING"].get(
                                    "experiments transferred", 0.0),
                            database=configuration[
                                "protein-protein interactions"]["STRING"].get(
                                    "database", 0.0),
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
                                    "combined score", 0.7),
                            taxon_identifier=configuration[
                                "protein-protein interactions"]["STRING"].get(
                                    "taxon identifier", 9606),
                            physical=configuration[
                                "protein-protein interactions"]["STRING"].get(
                                    "physical", False),
                        )

                    network.annotate_proteins()
                    network.remove_unannotated_proteins()

                    k += 1

                if "BioGRID" in configuration["protein-protein interactions"]:
                    biogrid.add_interactions(
                        network,
                        experimental_system=configuration[
                            "protein-protein interactions"]["BioGRID"].get(
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
                        experimental_system_type=configuration[
                            "protein-protein interactions"]["BioGRID"].get(
                                "experimental system type",
                                [
                                    "physical",
                                ],
                            ),
                        taxon_identifier=configuration[
                            "protein-protein interactions"]["BioGRID"].get(
                                "taxon identifier", 9606),
                        multi_validated_physical=configuration[
                            "protein-protein interactions"]["BioGRID"].get(
                                "multi-validated physical", False),
                    )

                if "CORUM" in configuration["protein-protein interactions"]:
                    corum.add_interactions(
                        network,
                        protein_complex_purification_method=configuration[
                            "protein-protein interactions"]["CORUM"].get(
                                "protein complex purification method",
                                [],
                            ),
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
                        ["IntAct"].get("MI score", 0.27),
                    )

                if "Reactome" in configuration["protein-protein interactions"]:
                    reactome.add_interactions(
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
                    string.add_interactions(
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
                                "experiments", 0.7),
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
                                "combined score", 0.7),
                        taxon_identifier=configuration[
                            "protein-protein interactions"]["STRING"].get(
                                "taxon identifier", 9606),
                        physical=configuration["protein-protein interactions"]
                        ["STRING"].get("physical", False),
                    )

            if "Cytoscape" in configuration:
                if (configuration["Cytoscape"].get(
                        "bar chart", {}).get("type") == "z-score"):
                    styles = CytoscapeStyles(
                        network,
                        bar_chart_range=configuration["Cytoscape"]
                        ["bar chart"].get("range", (-2.0, 2.0)),
                        get_bar_chart_range=network.get_z_score_range,
                        combine_sites=combine.COMBINE_SITES.get(
                            configuration["Cytoscape"]["bar chart"].get(
                                "combine sites"),
                            combine.COMBINE_SITES["absmax"],
                        ),
                    )

                elif (configuration["Cytoscape"].get(
                        "bar chart", {}).get("type") == "proportion"):
                    styles = CytoscapeStyles(
                        network,
                        bar_chart_range=configuration["Cytoscape"]
                        ["bar chart"].get("range", (0.025, 0.975)),
                        get_bar_chart_range=network.get_propotion_range,
                        combine_sites=combine.COMBINE_SITES.get(
                            configuration["Cytoscape"]["bar chart"].get(
                                "combine sites"),
                            combine.COMBINE_SITES["absmax"],
                        ),
                    )

                else:
                    styles = CytoscapeStyles(
                        network,
                        bar_chart_range=configuration["Cytoscape"]
                        ["bar chart"].get("range", (-1.0, 1.0)),
                        combine_sites=combine.COMBINE_SITES.get(
                            configuration["Cytoscape"]["bar chart"].get(
                                "combine sites"),
                            combine.COMBINE_SITES["absmax"],
                        ),
                    )

                export.export_styles(
                    styles,
                    os.path.splitext(os.path.basename(configuration_file))[0],
                    ".{}".format(i) if len(configurations) > 1 else "",
                )

                network.set_post_translational_modification()

                if (configuration["Cytoscape"].get(
                        "node color", {}).get("type") == "z-score"):
                    network.set_changes(
                        combine_sites=combine.COMBINE_SITES[
                            configuration["Cytoscape"]["node color"].get(
                                "combine sites", "absmax")],
                        changes=configuration["Cytoscape"]["node color"].get(
                            "range", (-2.0, 2.0)),
                        get_range=network.get_z_score_range,
                    )

                elif (configuration["Cytoscape"].get(
                        "node color", {}).get("type") == "proportion"):
                    network.set_changes(
                        combine_sites=combine.COMBINE_SITES[
                            configuration["Cytoscape"]["node color"].get(
                                "combine sites", "absmax")],
                        changes=configuration["Cytoscape"]["node color"].get(
                            "range", (0.025, 0.975)),
                        get_range=network.get_proportion_range,
                    )

                else:
                    network.set_changes(
                        combine_sites=combine.COMBINE_SITES[
                            configuration["Cytoscape"]["node color"].get(
                                "combine sites", "absmax")],
                        changes=configuration["Cytoscape"]["node color"].get(
                            "range", (-1.0, 1.0)),
                    )

            export.export_network(
                network,
                os.path.splitext(os.path.basename(configuration_file))[0],
                ".{}".format(i) if len(configurations) > 1 else "",
            )

            if "post-processing" in configuration:
                if "neighborhood extraction" in configuration:
                    for protein in configuration["neighborhood"].get(
                            "proteins", []):
                        export.export_network(
                            network.get_neighborhood(
                                protein,
                                configuration["neighborhood"].get(
                                    "distance", 1),
                                configuration["neighborhood"].get(
                                    "isoforms", True),
                            ),
                            os.path.splitext(
                                os.path.basename(configuration_file))[0],
                            ".{}.{}.{}".format(
                                i, protein,
                                network.nodes[protein]["gene name"])
                            if len(configurations) > 1 else ".{}.{}".format(
                                protein, network.nodes[protein]["gene name"]),
                        )

                if "module detection" in configuration["post-processing"]:
                    network.set_edge_weights(
                        weight=combine.COMBINE_CONFIDENCE_SCORES[
                            configuration["post-processing"]
                            ["module detection"].get("edge weight", "number")])

                    for j, module in enumerate(
                            network.get_modules(
                                module_size=configuration["post-processing"]
                                ["module detection"].get("module size", 0),
                                combine_module_sizes=combine.
                                COMBINE_MODULE_SIZES.get(
                                    configuration["post-processing"]
                                    ["module detection"].get("combine sizes"),
                                    combine.COMBINE_MODULE_SIZES["mean"],
                                ),
                                algorithm=algorithm.ALGORITHM.get(
                                    configuration["post-processing"]
                                    ["module detection"].get("algorithm"),
                                    algorithm.ALGORITHM["Louvain"],
                                ),
                            ),
                            start=1,
                    ):
                        export.export_network(
                            network.subgraph(module),
                            os.path.splitext(
                                os.path.basename(configuration_file))[0],
                            ".{}.{}".format(i, j)
                            if len(configurations) > 1 else ".{}".format(j),
                        )

                if "enrichment analysis" in configuration["post-processing"]:
                    network.set_edge_weights(
                        weight=combine.COMBINE_CONFIDENCE_SCORES[
                            configuration["post-processing"]
                            ["enrichment analysis"].get(
                                "edge weight", "number")])

                    if (configuration["post-processing"]
                        ["enrichment analysis"].get("type") == "z-score"):
                        modules, p_values = network.get_module_change_enrichment(
                            p=configuration["post-processing"]
                            ["enrichment analysis"].get("p", 0.05),
                            changes=configuration["post-processing"]
                            ["enrichment analysis"].get("range", (-2.0, 2.0)),
                            get_range=network.get_z_score_range,
                            combine_sites=combine.COMBINE_SITES.get(
                                configuration["post-processing"]
                                ["enrichment analysis"].get("combine sites"),
                                combine.COMBINE_SITES["absmax"],
                            ),
                            module_size=configuration["post-processing"]
                            ["enrichment analysis"].get("module size", 0),
                            combine_module_sizes=combine.COMBINE_MODULE_SIZES.
                            get(
                                configuration["post-processing"]
                                ["enrichment analysis"].get("combine sites"),
                                combine.COMBINE_MODULE_SIZES["mean"],
                            ),
                            test=configuration["post-processing"]
                            ["enrichment analysis"].get("test", "two-sided"),
                            algorithm=algorithm.ALGORITHM.get(
                                configuration["post-processing"]
                                ["enrichment analysis"].get("algorithm"),
                                algorithm.ALGORITHM["Louvain"],
                            ),
                        )

                    elif (configuration["post-processing"]
                          ["enrichment analysis"].get("type") == "proportion"):
                        modules, p_values = network.get_module_change_enrichment(
                            p=configuration["post-processing"]
                            ["enrichment analysis"].get("p", 0.05),
                            changes=configuration["post-processing"]
                            ["enrichment analysis"].get(
                                "range", (0.025, 0.975)),
                            get_range=network.get_proportion_range,
                            combine_sites=combine.COMBINE_SITES.get(
                                configuration["post-processing"]
                                ["enrichment analysis"].get("combine sites"),
                                combine.COMBINE_SITES["absmax"],
                            ),
                            module_size=configuration["post-processing"]
                            ["enrichment analysis"].get("module size", 0),
                            combine_module_sizes=combine.COMBINE_MODULE_SIZES.
                            get(
                                configuration["post-processing"]
                                ["enrichment analysis"].get("combine sites"),
                                combine.COMBINE_MODULE_SIZES["mean"],
                            ),
                            test=configuration["post-processing"]
                            ["enrichment analysis"].get("test", "two-sided"),
                            algorithm=algorithm.ALGORITHM.get(
                                configuration["post-processing"]
                                ["enrichment analysis"].get("algorithm"),
                                algorithm.ALGORITHM["Louvain"],
                            ),
                        )

                    else:
                        modules, p_values = network.get_module_change_enrichment(
                            p=configuration["post-processing"]
                            ["enrichment analysis"].get("p", 0.05),
                            changes=configuration["post-processing"]
                            ["enrichment analysis"].get("range", (-1.0, 1.0)),
                            combine_sites=combine.COMBINE_SITES.get(
                                configuration["post-processing"]
                                ["enrichment analysis"].get("combine sites"),
                                combine.COMBINE_SITES["absmax"],
                            ),
                            module_size=configuration["post-processing"]
                            ["enrichment analysis"].get("module size", 0),
                            combine_module_sizes=combine.COMBINE_MODULE_SIZES.
                            get(
                                configuration["post-processing"]
                                ["enrichment analysis"].get("combine sites"),
                                combine.COMBINE_MODULE_SIZES["mean"],
                            ),
                            test=configuration["post-processing"]
                            ["enrichment analysis"].get("test", "two-sided"),
                            algorithm=algorithm.ALGORITHM.get(
                                configuration["post-processing"]
                                ["enrichment analysis"].get("algorithm"),
                                algorithm.ALGORITHM["Louvain"],
                            ),
                        )

                    for time in p_values:
                        for ptm in sorted(p_values[time]):
                            for j in sorted(p_values[time][ptm]):
                                logger.info("{}\t{}\t{}\t{:.2E}".format(
                                    time,
                                    ptm,
                                    j + 1,
                                    p_values[time][ptm][j],
                                ))
                                export.export_network(
                                    network.subgraph(modules[j]),
                                    os.path.splitext(
                                        os.path.basename(
                                            configuration_file))[0],
                                    ".{}.{}".format(i, j + 1)
                                    if len(configurations) > 1 else
                                    ".{}".format(j + 1),
                                )


if __name__ == "__main__":
    main()
