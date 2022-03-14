import bisect
import math
import os
import re
import statistics

import networkx as nx
import pandas as pd

from enrichment import correction, test
from databases import biogrid, intact, mint, reactome, string, uniprot
from modularization import modularization


def get_protein_protein_interaction_network():
    return nx.Graph()


def annotate_proteins(network, taxon_identifier=9606):
    for accessions, gene_name, protein_name in uniprot.get_swissprot_entries(
            taxon_identifier):
        for protein in network:
            if protein.split("-")[0] == accessions[0]:
                network.nodes[protein]["gene"] = gene_name
                network.nodes[protein]["protein"] = protein_name


def remove_unannotated_proteins(network):
    network.remove_nodes_from([
        node for node in network if not (network.nodes[node].get("gene")
                                         or network.nodes[node].get("protein"))
    ])


def add_genes_from(network, genes, taxon_identifier=9606):
    for accessions, gene_name, protein_name in uniprot.get_swissprot_entries(
            taxon_identifier):
        if gene_name in genes:
            network.add_node(accessions[0])
            network.nodes[accessions[0]]["gene"] = gene_name
            network.nodes[accessions[0]]["protein"] = protein_name


def add_genes_from_table(
    network,
    file_name,
    gene_accession_column,
    gene_accession_format=re.compile("^(.+?)$"),
    sheet_name=0,
    header=0,
    taxon_identifier=9606,
):
    if os.path.splitext(file_name)[1].lstrip(".") in (
            "xls",
            "xlsx",
            "xlsm",
            "xlsb",
            "odf",
            "ods",
            "odt",
    ):
        table = pd.read_excel(
            file_name,
            sheet_name=sheet_name,
            header=header,
            usecols=[gene_accession_column],
            dtype={gene_accession_column: str},
        )
    else:
        table = pd.read_table(
            file_name,
            header=header,
            usecols=[gene_accession_column],
            dtype={gene_accession_column: str},
        )

    genes = set()
    for _, row in table.iterrows():
        if pd.isna(row[gene_accession_column]):
            continue

        genes.update([
            str(gene_accession) for gene_accession in
            gene_accession_format.findall(row[gene_accession_column])
        ])

    add_genes_from(network, genes, taxon_identifier)


def add_proteins_from(network, proteins):
    proteins_isoform = {}

    for protein_accession in proteins:
        if "-" in protein_accession and protein_accession.split(
                "-")[1].isnumeric():
            protein, isoform = protein_accession.split("-")
        else:
            protein, isoform = protein_accession, "0"

        if protein not in proteins_isoform:
            proteins_isoform[protein] = set()

        proteins_isoform[protein].add(isoform)

    primary_accession = {}
    for accessions, gene_name, protein_name in uniprot.get_swissprot_entries():
        for i, accession in enumerate(accessions):
            if accession in proteins_isoform:
                for isoform in proteins_isoform[accession]:
                    if isoform == "0":
                        network.add_node(accessions[0])
                        network.nodes[accessions[0]]["gene"] = gene_name
                        network.nodes[accessions[0]]["protein"] = protein_name

                    elif i == 0:
                        network.add_node("{}-{}".format(accession, isoform))
                        network.nodes["{}-{}".format(
                            accession, isoform)]["gene"] = gene_name
                        network.nodes["{}-{}".format(
                            accession, isoform)]["protein"] = protein_name

                if i > 0:
                    if accession not in primary_accession:
                        primary_accession[accession] = set()

                    primary_accession[accession].add(accessions[0])

    return primary_accession


def add_proteins_from_table(network,
                            file_name,
                            protein_accession_column,
                            protein_accession_format=re.compile("^(.+?)$"),
                            time=0,
                            modification="",
                            position_column="",
                            position_format=re.compile("^(.+?)$"),
                            replicates=[],
                            sheet_name=0,
                            header=0,
                            num_sites=100,
                            num_replicates=1,
                            replicate_combination=statistics.mean,
                            measurement_conversion=math.log2):
    if os.path.splitext(file_name)[1].lstrip(".") in (
            "xls",
            "xlsx",
            "xlsm",
            "xlsb",
            "odf",
            "ods",
            "odt",
    ):
        table = pd.read_excel(
            file_name,
            sheet_name=sheet_name,
            header=header,
            usecols=[protein_accession_column] +
            [column for column in (position_column, ) if column] + replicates,
            dtype={
                protein_accession_column: str,
                **{column: str
                   for column in (position_column, ) if column},
                **{column: float
                   for column in replicates},
            },
        )
    else:
        table = pd.read_table(
            file_name,
            header=header,
            usecols=[protein_accession_column] +
            [column for column in (position_column, ) if column] + replicates,
            dtype={
                protein_accession_column: str,
                **{column: str
                   for column in (position_column, ) if column},
                **{column: float
                   for column in replicates},
            },
        )

    proteins = {}
    for _, row in table.iterrows():
        if pd.isna(row[protein_accession_column]):
            continue

        protein_accessions = [
            str(protein_accession) for protein_accession in
            protein_accession_format.findall(row[protein_accession_column])
        ]

        if replicates:
            if position_column and not pd.isna(row[position_column]):
                positions = [
                    int(position) for position in position_format.findall(
                        row[position_column])
                ]
            else:
                positions = []

            if len(protein_accessions) != len(positions):
                if len(protein_accessions) > len(positions):
                    positions.extend([
                        0
                        for _ in range(len(positions), len(protein_accessions))
                    ])
                else:
                    positions = positions[:len(protein_accessions)]

            measurements = [
                row[replicate] for replicate in replicates
                if not pd.isna(row[replicate])
            ]

            if len(measurements) >= min(num_replicates, len(replicates)):
                for protein_accession, position in zip(protein_accessions,
                                                       positions):
                    if protein_accession not in proteins:
                        proteins[protein_accession] = []

                    bisect.insort(
                        proteins[protein_accession],
                        (
                            position,
                            measurement_conversion(
                                replicate_combination(measurements)),
                        ),
                    )
        else:
            for protein_accession in protein_accessions:
                if protein_accession not in proteins:
                    proteins[protein_accession] = []

    primary_accession = add_proteins_from(network, proteins.keys())

    for protein in list(proteins):
        if protein in primary_accession:
            for accession in sorted(primary_accession[protein]):
                if accession not in proteins:
                    proteins[accession] = proteins[protein]
            del proteins[protein]

    for protein in proteins:
        if protein not in network:
            continue

        proteins[protein] = sorted(
            sorted(
                proteins[protein],
                key=lambda item: abs(item[1]),
            )[-num_sites:])

        for i in range(len(proteins[protein])):
            network.nodes[protein]["{}-{}-{}".format(
                time, modification, i + 1)] = proteins[protein][i][1]


def get_times(network):
    return tuple(
        sorted(
            set(
                int(change.split("-")[0]) for protein in network
                for change in network.nodes[protein]
                if len(change.split("-")) == 3
                and change.split("-")[0].isnumeric())))


def get_post_translational_modifications(network, time):
    return tuple(
        sorted(
            set(
                change.split("-")[1] for protein in network
                for change in network.nodes[protein]
                if len(change.split("-")) == 3
                and change.split("-")[0] == str(time))))


def get_sites(network, time, modification):
    return max(
        int(change.split("-")[2]) for protein in network
        for change in network.nodes[protein]
        if len(change.split("-")) == 3 and change.split("-")[0] == str(time)
        and change.split("-")[1] == modification)


def get_changes(network,
                time,
                modification,
                site_combination=lambda sites: max(sites, key=abs)):
    changes = []
    for protein in network:
        sites = [
            network.nodes[protein][change] for change in network.nodes[protein]
            if len(change.split("-")) == 3 and change.split("-")[0] == str(
                time) and change.split("-")[1] == modification
        ]

        if sites:
            changes.append(site_combination(sites))

    return changes


def get_standard_score_of_change(
        network,
        change,
        time,
        modification,
        site_combination=lambda sites: max(sites, key=abs),
):
    changes = get_changes(network, time, modification, site_combination)
    mean = statistics.mean(changes)
    stdev = statistics.stdev(changes, xbar=mean)
    return (change - mean) / stdev


def get_change_of_standard_score(
    network,
    standard_score,
    time,
    modification,
    site_combination=lambda sites: max(sites, key=abs)):
    changes = get_changes(network, time, modification, site_combination)
    mean = statistics.mean(changes)
    stdev = statistics.stdev(changes, xbar=mean)
    return standard_score * stdev + mean


def get_quantile_of_change(
        network,
        change,
        time,
        modification,
        site_combination=lambda sites: max(sites, key=abs),
):
    changes = get_changes(network, time, modification, site_combination)
    return len([c for c in changes if c <= change]) / len(changes)


def get_change_of_quantile(network,
                           quantile,
                           time,
                           modification,
                           site_combination=lambda sites: max(sites, key=abs)):
    changes = get_changes(network, time, modification, site_combination)
    return sorted(changes)[math.floor(quantile * (len(changes) - 1))]


def set_post_translational_modification(network):
    for time in get_times(network):
        for protein in network:
            network.nodes[protein]["post-translational modification {}".format(
                time)] = "".join(
                    sorted(
                        set(
                            change.split("-")[1]
                            for change in network.nodes[protein]
                            if len(change.split("-")) == 3
                            and change.split("-")[0] == str(time))))


def set_changes(
    network,
    site_combination=lambda sites: max(sites, key=abs),
    changes=(-1.0, 1.0),
    get_change=lambda network, change, time, modification, site_combination:
    change,
):
    times = get_times(network)
    modifications = {
        time: get_post_translational_modifications(network, time)
        for time in times
    }

    change_range = {
        time: {
            modification: (get_change(
                network,
                changes[0],
                time,
                modification,
                site_combination,
            ),
                           get_change(
                               network,
                               changes[1],
                               time,
                               modification,
                               site_combination,
                           ))
            for modification in modifications[time]
        }
        for time in times
    }

    for time in times:
        for protein in network:
            classification = {}
            for modification in modifications[time]:
                sites = [
                    network.nodes[protein][change]
                    for change in network.nodes[protein]
                    if len(change.split("-")) == 3 and change.split("-")[0] ==
                    str(time) and change.split("-")[1] == modification
                ]

                if sites:
                    combined_change = site_combination(sites)

                    if combined_change >= 1.0 * change_range[time][
                            modification][1]:
                        classification[modification] = "up"
                    elif 1.0 * change_range[time][modification][
                            1] > combined_change >= 0.5 * change_range[time][
                                modification][1]:
                        classification[modification] = "mid_up"
                    elif 0.5 * change_range[time][modification][
                            1] > combined_change > 0.5 * change_range[time][
                                modification][0]:
                        classification[modification] = "mid"
                    elif 0.5 * change_range[time][modification][
                            0] >= combined_change > 1.0 * change_range[time][
                                modification][0]:
                        classification[modification] = "mid_down"
                    else:
                        classification[modification] = "down"

            if classification:
                if set(classification.values()) == {"up"}:
                    network.nodes[protein]["change {}".format(time)] = "up"
                elif set(classification.values()) == {"mid_up", "up"}:
                    network.nodes[protein]["change {}".format(time)] = "up"
                elif set(classification.values()) == {"mid", "mid_up", "up"}:
                    network.nodes[protein]["change {}".format(time)] = "up"

                elif set(classification.values()) == {"mid_up"}:
                    network.nodes[protein]["change {}".format(time)] = "mid_up"
                elif set(classification.values()) == {"mid", "mid_up"}:
                    network.nodes[protein]["change {}".format(time)] = "mid_up"

                elif set(classification.values()) == {"mid"}:
                    network.nodes[protein]["change {}".format(time)] = "mid"

                elif set(classification.values()) == {"mid", "mid_down"}:
                    network.nodes[protein]["change {}".format(
                        time)] = "mid_down"
                elif set(classification.values()) == {"mid_down"}:
                    network.nodes[protein]["change {}".format(
                        time)] = "mid_down"

                elif set(classification.values()) == {
                        "mid", "mid_down", "down"
                }:
                    network.nodes[protein]["change {}".format(time)] = "down"
                elif set(classification.values()) == {"mid_down", "down"}:
                    network.nodes[protein]["change {}".format(time)] = "down"
                elif set(classification.values()) == {"down"}:
                    network.nodes[protein]["change {}".format(time)] = "down"
                else:
                    network.nodes[protein]["change {}".format(
                        time)] = " ".join(
                            "{}_{}".format(modification,
                                           classification[modification])
                            for modification in sorted(classification.keys()))

            else:
                network.nodes[protein]["change {}".format(time)] = "mid"


def add_proteins_from_biogrid(network,
                              experimental_system=[],
                              experimental_system_type=[],
                              interaction_throughput=[],
                              multi_validated_physical=False,
                              taxon_identifier=9606):
    nodes_to_add = set()
    for interactor_a, interactor_b in biogrid.get_proteins(
            experimental_system, experimental_system_type,
            interaction_throughput, multi_validated_physical,
            taxon_identifier):
        if (interactor_a in network and interactor_b not in network):
            nodes_to_add.add(interactor_b)

        elif (interactor_a not in network and interactor_b in network):
            nodes_to_add.add(interactor_a)

    network.add_nodes_from(nodes_to_add)


def add_protein_protein_interactions_from_biogrid(
    network,
    experimental_system=[],
    experimental_system_type=[],
    interaction_throughput=[],
    multi_validated_physical=False,
    taxon_identifier=9606,
):
    for interactor_a, interactor_b in biogrid.get_protein_protein_interactions(
            experimental_system, experimental_system_type,
            interaction_throughput, multi_validated_physical,
            taxon_identifier):
        if (interactor_a in network and interactor_b in network
                and interactor_a != interactor_b):
            network.add_edge(interactor_a, interactor_b)
            network.edges[interactor_a, interactor_b]["BioGRID"] = 1.0


def add_proteins_from_intact(network,
                             interaction_detection_methods=[],
                             interaction_types=[],
                             mi_score=0.0,
                             taxon_identifier=9606):
    nodes_to_add = set()
    for interactor_a, interactor_b in intact.get_proteins(
            interaction_detection_methods, interaction_types, mi_score,
            taxon_identifier):
        if (interactor_a in network and interactor_b not in network):
            nodes_to_add.add(interactor_b)

        elif (interactor_a not in network and interactor_b in network):
            nodes_to_add.add(interactor_a)

    network.add_nodes_from(nodes_to_add)


def add_protein_protein_interactions_from_intact(
        network,
        interaction_detection_methods=[],
        interaction_types=[],
        mi_score=0.0,
        taxon_identifier=9606):

    for interactor_a, interactor_b, score in intact.get_protein_protein_interactions(
            interaction_detection_methods, interaction_types, mi_score,
            taxon_identifier):
        if (interactor_a in network and interactor_b in network
                and interactor_a != interactor_b):
            if network.has_edge(interactor_a, interactor_b):
                network.edges[interactor_a, interactor_b]["IntAct"] = max(
                    score, network.edges[interactor_a,
                                         interactor_b].get("IntAct", 0.0))
            else:
                network.add_edge(interactor_a, interactor_b)
                network.edges[interactor_a, interactor_b]["IntAct"] = score


def add_proteins_from_mint(network,
                           interaction_detection_methods=[],
                           interaction_types=[],
                           mi_score=0.0,
                           taxon_identifier=9606):

    nodes_to_add = set()
    for interactor_a, interactor_b in mint.get_proteins(
            interaction_detection_methods, interaction_types, mi_score,
            taxon_identifier):
        if (interactor_a in network and interactor_b not in network):
            nodes_to_add.add(interactor_b)

        elif (interactor_a not in network and interactor_b in network):
            nodes_to_add.add(interactor_a)

    network.add_nodes_from(nodes_to_add)


def add_protein_protein_interactions_from_mint(
        network,
        interaction_detection_methods=[],
        interaction_types=[],
        mi_score=0.0,
        taxon_identifier=9606):

    for interactor_a, interactor_b, score in mint.get_protein_protein_interactions(
            interaction_detection_methods, interaction_types, mi_score,
            taxon_identifier):
        if (interactor_a in network and interactor_b in network
                and interactor_a != interactor_b):
            if network.has_edge(interactor_a, interactor_b):
                network.edges[interactor_a, interactor_b]["MINT"] = max(
                    score, network.edges[interactor_a,
                                         interactor_b].get("MINT", 0.0))
            else:
                network.add_edge(interactor_a, interactor_b)
                network.edges[interactor_a, interactor_b]["MINT"] = score


def add_proteins_from_reactome(network,
                               interaction_type=[],
                               interaction_context=[],
                               taxon_identifier=9606):
    nodes_to_add = set()

    for interactor_a, interactor_b in reactome.get_proteins(
            interaction_type, interaction_context, taxon_identifier):
        if (interactor_a in network and interactor_b not in network):
            nodes_to_add.add(interactor_b)

        elif (interactor_a not in network and interactor_b in network):
            nodes_to_add.add(interactor_a)

    network.add_nodes_from(nodes_to_add)


def add_protein_protein_interactions_from_reactome(
    network,
    interaction_type=[],
    interaction_context=[],
    taxon_identifier=9606,
):
    for interactor_a, interactor_b in reactome.get_protein_protein_interactions(
            interaction_type, interaction_context, taxon_identifier):
        if (interactor_a in network and interactor_b in network
                and interactor_a != interactor_b):
            network.add_edge(interactor_a, interactor_b)
            network.edges[interactor_a, interactor_b]["Reactome"] = 1.0


def add_proteins_from_string(network,
                             neighborhood=0.0,
                             neighborhood_transferred=0.0,
                             fusion=0.0,
                             cooccurence=0.0,
                             homology=0.0,
                             coexpression=0.0,
                             coexpression_transferred=0.0,
                             experiments=0.0,
                             experiments_transferred=0.0,
                             database=0.0,
                             database_transferred=0.0,
                             textmining=0.0,
                             textmining_transferred=0.0,
                             combined_score=0.0,
                             physical=False,
                             taxon_identifier=9606,
                             version=11.5):
    nodes_to_add = set()

    for interactor_a, interactor_b in string.get_proteins(
            neighborhood, neighborhood_transferred, fusion, cooccurence,
            homology, coexpression, coexpression_transferred, experiments,
            experiments_transferred, database, database_transferred,
            textmining, textmining_transferred, combined_score, physical,
            taxon_identifier, version):
        if (interactor_a in network and interactor_b not in network):
            nodes_to_add.add(interactor_b)

        elif (interactor_a not in network and interactor_b in network):
            nodes_to_add.add(interactor_a)

    network.add_nodes_from(nodes_to_add)


def add_protein_protein_interactions_from_string(network,
                                                 neighborhood=0.0,
                                                 neighborhood_transferred=0.0,
                                                 fusion=0.0,
                                                 cooccurence=0.0,
                                                 homology=0.0,
                                                 coexpression=0.0,
                                                 coexpression_transferred=0.0,
                                                 experiments=0.0,
                                                 experiments_transferred=0.0,
                                                 database=0.0,
                                                 database_transferred=0.0,
                                                 textmining=0.0,
                                                 textmining_transferred=0.0,
                                                 combined_score=0.0,
                                                 physical=False,
                                                 taxon_identifier=9606,
                                                 version=11.5):
    for interactor_a, interactor_b, score in string.get_protein_protein_interactions(
            neighborhood, neighborhood_transferred, fusion, cooccurence,
            homology, coexpression, coexpression_transferred, experiments,
            experiments_transferred, database, database_transferred,
            textmining, textmining_transferred, combined_score, physical,
            taxon_identifier, version):
        if (interactor_a in network and interactor_b in network
                and interactor_a != interactor_b):
            if network.has_edge(interactor_a, interactor_b):
                network.edges[interactor_a, interactor_b]["STRING"] = max(
                    score,
                    network.edges[interactor_a,
                                  interactor_b].get("STRING", 0.0),
                )
            else:
                network.add_edge(interactor_a, interactor_b)
                network.edges[interactor_a, interactor_b]["STRING"] = (score)


def get_databases(network):
    return tuple(
        sorted(
            set(database for edge in network.edges()
                for database in network.edges[edge]).intersection(
                    {"BioGRID", "IntAct", "MINT", "Reactome", "STRING"})))


def set_edge_weights(
    network,
    weight=lambda confidence_scores: int(bool(confidence_scores)),
    attribute="weight",
):
    databases = get_databases(network)

    for edge in network.edges:
        network.edges[edge][attribute] = weight({
            database: network.edges[edge][database]
            for database in databases if database in network.edges[edge]
        })


def remove_edge_weights(network, attribute="weight"):
    for _, _, data in network.edges(data=True):
        if attribute in data:
            del data[attribute]


def get_modules(network,
                module_size,
                module_size_combination=statistics.mean,
                algorithm=modularization.louvain,
                resolution=1.0,
                weight="weight"):
    G = network.copy()
    G.remove_nodes_from(list(nx.isolates(G)))

    communities = algorithm(G, resolution, weight)

    while (module_size_combination(
            len(community) for community in communities) > module_size):

        subdivision = False
        for i, subdivided_community in enumerate(
                algorithm(G.subgraph(communities[j]), resolution, weight)
                for j in range(len(communities))):
            if len(subdivided_community) > 1:
                subdivision = True
                communities[i:i + 1] = subdivided_community

        if not subdivision:
            break

    return [network.subgraph(community) for community in communities]


def get_proteins(
    network,
    time,
    modification,
    site_combination=lambda sites: max(sites, key=abs),
    combined_change_filter=lambda combined_sites: bool(combined_sites)):
    proteins = []
    for protein in network:
        sites = [
            network.nodes[protein][change] for change in network.nodes[protein]
            if len(change.split("-")) == 3 and change.split("-")[0] == str(
                time) and change.split("-")[1] == modification
        ]

        if sites and combined_change_filter(site_combination(sites)):
            proteins.append(protein)

    return proteins


def get_binary_change_enrichment(
    network,
    modules,
    changes=(-1.0, 1.0),
    get_change=lambda network, change, time, modification, site_combination:
    change,
    site_combination=lambda sites: max(sites, key=abs),
    test=test.hypergeometric,
    correction=correction.benjamini_hochberg,
):
    p_values = {}
    for time in get_times(network):
        for modification in get_post_translational_modifications(
                network, time):
            modified_proteins = len(get_proteins(network, time, modification))

            modified_module_proteins = [
                len(get_proteins(module, time, modification))
                for module in modules
            ]

            change_range = (get_change(network, changes[0], time, modification,
                                       site_combination),
                            get_change(network, changes[1], time, modification,
                                       site_combination))

            target_proteins = len(
                get_proteins(
                    network,
                    time,
                    modification,
                    site_combination=site_combination,
                    combined_change_filter=lambda combined_change:
                    combined_change <= change_range[
                        0] or combined_change >= change_range[1],
                ))

            target_module_proteins = [
                len(
                    get_proteins(
                        module,
                        time,
                        modification,
                        site_combination=site_combination,
                        combined_change_filter=lambda combined_change:
                        combined_change <= change_range[
                            0] or combined_change >= change_range[1],
                    )) for module in modules
            ]

            p_values.update({(module, time, modification):
                             test(target_module_proteins[i], modified_proteins,
                                  target_proteins, modified_module_proteins[i])
                             for i, module in enumerate(modules)})

    p_values = correction(p_values)

    return {
        module: {
            time: {
                modification: p_values[(module, time, modification)]
                for modification in get_post_translational_modifications(
                    network, time)
            }
            for time in get_times(network)
        }
        for module in modules
    }


def get_continuos_change_enrichment(
    network,
    modules,
    site_combination=lambda sites: max(sites, key=abs),
    test=test.wilcoxon,
    correction=correction.benjamini_hochberg,
):
    p_values = {}
    for time in get_times(network):
        for modification in get_post_translational_modifications(
                network, time):

            network_changes = [
                get_changes(nx.union_all([m for m in modules if m != module]),
                            time, modification, site_combination)
                for module in modules
            ]

            module_changes = [
                get_changes(module, time, modification, site_combination)
                for module in modules
            ]

            p_values.update({(module, time, modification):
                             test(module_changes[i], network_changes[i])
                             for i, module in enumerate(modules)})

    p_values = correction(p_values)

    return {
        module: {
            time: {
                modification: p_values[(module, time, modification)]
                for modification in get_post_translational_modifications(
                    network, time)
            }
            for time in get_times(network)
        }
        for module in modules
    }


def export(network, basename, suffix=""):
    nx.write_graphml_xml(network,
                         "{0}{1}.graphml".format(basename, suffix),
                         named_key_ids=True,
                         infer_numeric_types=True)
