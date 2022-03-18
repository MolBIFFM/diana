"""protein-protein interaction network"""
import bisect
import math
import os
import re
import statistics
from typing import Callable, Collection, Container, Hashable, Optional, Union

import networkx as nx
import pandas as pd
from analysis import correction, modularization, test
from databases import (biogrid, gene_ontology, intact, mint, reactome, string,
                       uniprot)


def get_protein_protein_interaction_network() -> nx.Graph():
    """
    Initializes a protein-protein interaction network.

    Returns
        The protein-protein interaction network.

    """
    return nx.Graph()


def annotate_proteins(network: nx.Graph,
                      taxonomy_identifier: int = 9606) -> None:
    """
    Associate nodes with gene and protein names from SwissProt.

    Args:
        network: The protein-protein interaction network.
        taxonomy_identifier: The taxonomy identifier.
    """
    for accessions, gene_name, protein_name in uniprot.get_swissprot_entries(
            taxonomy_identifier):
        for protein in network:
            if protein.split("-")[0] == accessions[0]:
                network.nodes[protein]["gene"] = gene_name
                network.nodes[protein]["protein"] = protein_name


def remove_unannotated_proteins(network: nx.Graph) -> None:
    """
    Remove proteins not associated with a gene or protein name.

    Args:
        network: The protein-protein interaction network.
    """
    network.remove_nodes_from([
        node for node in network if not (network.nodes[node].get("gene") or
                                         network.nodes[node].get("protein"))
    ])


def add_genes_from(network: nx.Graph,
                   genes: Container,
                   taxonomy_identifier: int = 9606) -> None:
    """
    Add primary Uniprot protein accessions corresponding to genes to a
    protein-protein interaction network.

    Args:
        network: The protein-protein interaction network.
        genes: The Uniprot gene accessions whose corresponding primary Uniprot
            protein accessions are added.
        taxonomy_identifier: The taxonomy identifier.
    """
    for accessions, gene_name, protein_name in uniprot.get_swissprot_entries(
            taxonomy_identifier):
        if gene_name in genes:
            network.add_node(accessions[0])
            network.nodes[accessions[0]]["gene"] = gene_name
            network.nodes[accessions[0]]["protein"] = protein_name


def add_genes_from_table(
    network: nx.Graph,
    file_name: str,
    gene_accession_column: Union[int, str],
    gene_accession_format: re.Pattern = re.compile("^(.+?)$"),
    sheet_name=Union[int, str],
    header: int = 0,
    taxonomy_identifier: int = 9606,
) -> None:
    """
    Parse UniProt gene accessions from a tabular file and add the corresponding
    primary Uniprot protein accessions to a protein-protein interaction network.

    Args:
        network: The protein-protein interaction network.
        file_name: The file location of the file.
        gene_accession_column: The column containing UniProt gene accessions.
        gene_accession_format: A regular expression to extract gene accessions
            from an entry.
        sheet_name: The sheet to parse gene accessions from.
        header: The index of the header row.
        taxonomy_identifier: The taxonomy identifier.
    """
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

    add_genes_from(network, genes, taxonomy_identifier)


def add_proteins_from(network: nx.Graph, proteins: Container) -> None:
    """
    Add primary Uniprot protein accessions corresponding to proteins to a
    protein-protein interaction network.

    Args:
        network: The protein-protein interaction network.
        proteins: The protein accessions whose corresponding primary Uniprot
            protein accessions are added.
    """
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


def add_proteins_from_table(
    network: nx.Graph,
    file_name: str,
    protein_accession_column: Union[int, str],
    protein_accession_format: re.Pattern = re.compile("^(.+?)$"),
    time: int = 0,
    modification: str = "",
    position_column: str = Union[int, str],
    position_format: re.Pattern = re.compile("^(.+?)$"),
    replicates: Optional[Container[Union[int, str]]] = None,
    sheet_name: Union[int, str] = 0,
    header: int = 0,
    num_sites: int = 100,
    num_replicates: int = 1,
    replicate_combination: Callable[[Container[float]],
                                    float] = statistics.mean,
    measurement_conversion: Callable[
        [Container[float], Callable[[Container[float]], float]],
        float] = lambda changes, combination: math.log2(combination(changes))
) -> None:
    """
    Parse UniProt protein accessions and measurements of changes from a tabular
    file and add the corresponding primary Uniprot protein accessions to a
    protein-protein interaction network.

    Args:
        network: The protein-protein interaction network.
        file_name: The file location of the file.
        protein_accession_column: The column containing UniProt protein
            accessions.
        protein_accession_format: A regular expression to extract protein
            accessions from a corresponding entry.
        time: The time of measurement to associate with changes.
        modification: An identifier for the type of post-translational
            modification to associate with changes.
        position_column: The column containing sites corresponding to
            measurements.
        replicates: The columns containing replicates of measurements.
        sheet_name: The sheet to parse protein accessions from.
        header: The index of the header row.
        num_sites: The maximum number of measurements to associate with a
            protein-accession, prioritized by largest absolute value.
        num_replicates: The minimum number of replicates to accept a
            measurement.
        replicate_combination: A function to combine replicates into a single
            change.
        measurement_conversion: A function to convert the measurements reported
            to log2-scale.
    """
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
            [column for column in (position_column,) if column] + replicates,
            dtype={
                protein_accession_column: str,
                **{column: str for column in (position_column,) if column},
                **{column: float for column in replicates},
            },
        )
    else:
        table = pd.read_table(
            file_name,
            header=header,
            usecols=[protein_accession_column] +
            [column for column in (position_column,) if column] + replicates,
            dtype={
                protein_accession_column: str,
                **{column: str for column in (position_column,) if column},
                **{column: float for column in replicates},
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
                row[replicate]
                for replicate in replicates
                if not pd.isna(row[replicate])
            ]

            if len(measurements) >= min(num_replicates, len(replicates)):
                for protein_accession, position in zip(protein_accessions,
                                                       positions):
                    if protein_accession not in proteins:
                        proteins[protein_accession] = []

                    bisect.insort(
                        proteins[protein_accession],
                        (position,
                         measurement_conversion(measurements,
                                                replicate_combination)),
                    )
        else:
            for protein_accession in protein_accessions:
                if protein_accession not in proteins:
                    proteins[protein_accession] = []

    primary_accession = add_proteins_from(network, proteins)

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
            network.nodes[protein]["{} {} {}".format(
                time, modification, i + 1)] = proteins[protein][i][1]


def get_proteins(
        network: nx.Graph,
        time: int,
        modification: str,
        site_combination: Optional[Callable[[Collection[float]], float]] = None,
        combined_change_filter: Callable[[float], bool] = bool) -> set[str]:
    """
    Returns proteins of a protein-protein interaction network with specified
    changes.

    Args:
        network: The protein-protein interaction network.
        time: The time of measurement.
        modification: The type of post-translational modification.
        site_combination: An optional function to derive protein-specific
            changes from site-specific changes.
        combined_change_filter: The predicate to filter combined changes.

    Returns:
        Proteins whose combined change of a particular type of
        post-translational modification at a particular time of measurement
        fulfil the predicate.
    """
    proteins = []
    for protein in network:
        sites = [
            network.nodes[protein][change]
            for change in network.nodes[protein]
            if len(change.split(" ")) == 3 and change.split(" ")[0] == str(time)
            and change.split(" ")[1] == modification
        ]

        if sites and combined_change_filter(site_combination(sites)):
            proteins.append(protein)

    return set(proteins)


def get_times(network: nx.Graph) -> tuple[int, ...]:
    """
    Returns the times of measurement represented in a protein-protein
    interaction network.

    Args:
        network: The protein-protein interaction network.

    Returns:
        The times of measurement associated with any changes of any protein in
        the protein-protein interaction network.
    """
    return tuple(
        sorted(
            set(
                int(change.split(" ")[0])
                for protein in network
                for change in network.nodes[protein]
                if len(change.split(" ")) == 3 and
                change.split(" ")[0].isnumeric())))


def get_post_translational_modifications(network: nx.Graph,
                                         time: int) -> tuple[str, ...]:
    """
    Returns the types of post-translational modification represented in a
    protein-protein interaction network at a particular time of measurement.

    Args:
        network: The protein-protein interaction network.
        time: The time of measurement.

    Returns:
        The types of post-translational modification associated with any changes
        of any protein in the protein-protein interaction network at a
        particular time of measurement.
    """
    return tuple(
        sorted(
            set(
                change.split(" ")[1]
                for protein in network
                for change in network.nodes[protein]
                if len(change.split(" ")) == 3 and
                change.split(" ")[0] == str(time))))


def get_sites(network: nx.Graph, time: int, modification: str) -> int:
    """
    Returns the number of sites of proteins in a protein-protein interaction
    network for a particular type of post-translational modification at a
    particular time of measurement.

    Args:
        network: The protein-protein interaction network.
        time: The time of measurement.
        modification: The type of post-translational modification.

    Returns:
        The maximum number of sites of any protein in a protein-protein
        interaction network for a particular type of post-translational
        modification at a particular time of measurement.
    """
    return max(
        int(change.split(" ")[2])
        for protein in network
        for change in network.nodes[protein]
        if len(change.split(" ")) == 3 and change.split(" ")[0] == str(time) and
        change.split(" ")[1] == modification)


def set_post_translational_modification(network: nx.Graph) -> None:
    """
    Annotate proteins with summary of type of corresponding modifications.

    Args:
        network: The protein-protein interaction network.
    """
    for time in get_times(network):
        for protein in network:
            network.nodes[protein]["post-translational modification {}".format(
                time)] = "".join(
                    sorted(
                        set(
                            change.split(" ")[1]
                            for change in network.nodes[protein]
                            if len(change.split(" ")) == 3 and
                            change.split(" ")[0] == str(time))))


def set_changes(
    network: nx.Graph,
    site_combination: Optional[Callable[[Collection[float]], float]] = None,
    changes: tuple[float, float] = (-1.0, 1.0),
    convert_change: Callable[
        [nx.Graph, float, int, str, Callable[[Collection[float]],
                                             float]], float] = lambda network,
    change, time, modification, site_combination: change,
) -> None:
    """
    Annotate nodes with summary of change trend.

    Args:
        network: The protein-protein interaction network.
        site_combination: An optional function to derive protein-specific changes from
            site-specific changes.
        changes: Proteins are classified by whether their representative exceed
            either this range, the range defined by half the bounds or none.
        convert_change: The function used convert the bounds to log2-fold
            changes.
    """
    times = get_times(network)
    modifications = {
        time: get_post_translational_modifications(network, time)
        for time in times
    }

    change_range = {
        time: {
            modification: (convert_change(
                network,
                changes[0],
                time,
                modification,
                site_combination,
            ),
                           convert_change(
                               network,
                               changes[1],
                               time,
                               modification,
                               site_combination,
                           )) for modification in modifications[time]
        } for time in times
    }

    for time in times:
        for protein in network:
            classification = {}
            for modification in modifications[time]:
                sites = [
                    network.nodes[protein][change]
                    for change in network.nodes[protein]
                    if len(change.split(" ")) == 3 and change.split(" ")[0] ==
                    str(time) and change.split(" ")[1] == modification
                ]

                if sites:
                    combined_change = site_combination(sites)

                    if combined_change >= 1.0 * change_range[time][
                            modification][1]:
                        classification[modification] = "up"
                    elif 1.0 * change_range[time][modification][
                            1] > combined_change >= 0.5 * change_range[time][
                                modification][1]:
                        classification[modification] = "mid up"
                    elif 0.5 * change_range[time][modification][
                            1] > combined_change > 0.5 * change_range[time][
                                modification][0]:
                        classification[modification] = "mid"
                    elif 0.5 * change_range[time][modification][
                            0] >= combined_change > 1.0 * change_range[time][
                                modification][0]:
                        classification[modification] = "mid down"
                    else:
                        classification[modification] = "down"

            if classification:
                if set(classification.values()) == {"up"}:
                    network.nodes[protein]["change {}".format(time)] = "up"
                elif set(classification.values()) == {"mid up", "up"}:
                    network.nodes[protein]["change {}".format(time)] = "up"
                elif set(classification.values()) == {"mid", "mid up", "up"}:
                    network.nodes[protein]["change {}".format(time)] = "up"

                elif set(classification.values()) == {"mid up"}:
                    network.nodes[protein]["change {}".format(time)] = "mid up"
                elif set(classification.values()) == {"mid", "mid up"}:
                    network.nodes[protein]["change {}".format(time)] = "mid up"

                elif set(classification.values()) == {"mid"}:
                    network.nodes[protein]["change {}".format(time)] = "mid"

                elif set(classification.values()) == {"mid", "mid down"}:
                    network.nodes[protein]["change {}".format(
                        time)] = "mid down"
                elif set(classification.values()) == {"mid down"}:
                    network.nodes[protein]["change {}".format(
                        time)] = "mid down"

                elif set(
                        classification.values()) == {"mid", "mid down", "down"}:
                    network.nodes[protein]["change {}".format(time)] = "down"
                elif set(classification.values()) == {"mid down", "down"}:
                    network.nodes[protein]["change {}".format(time)] = "down"
                elif set(classification.values()) == {"down"}:
                    network.nodes[protein]["change {}".format(time)] = "down"
                else:
                    network.nodes[protein]["change {}".format(time)] = " ".join(
                        "{}_{}".format(modification,
                                       classification[modification])
                        for modification in sorted(classification))

            else:
                network.nodes[protein]["change {}".format(time)] = "mid"


def add_proteins_from_biogrid(
        network: nx.Graph,
        experimental_system: Optional[Container[str]] = None,
        experimental_system_type: Optional[Container[str]] = None,
        interaction_throughput: Optional[Container[str]] = None,
        multi_validated_physical: bool = False,
        taxonomy_identifier: int = 9606,
        version: Optional[tuple[int, int, int]] = None) -> None:
    """
    Add proteins interacting with proteins in a protein-protein interaction
    network from BioGRID to the network.

    Args:
        network: The protein-protein interaction network.
        experimental_system: The accepted experimental evidence codes. If none
            are specified, any is accepted.
        experimental_system_type: The accepted categories of experimental
            evidence. If none are specified, any is accepted.
        interaction_throughput:  The accepted levels of interaction throughput.
            If none are specified, any is accepted.
        multi-validated physical: If True, consider only multi-validated
            physical interactions.
        taxonomy_identifier: The taxonomy identifier.
        version: The version of the BioGRID database, if not passed, the latest.
    """
    nodes_to_add = set()
    for interactor_a, interactor_b in biogrid.get_protein_protein_interactions(
            experimental_system, experimental_system_type,
            interaction_throughput, multi_validated_physical,
            taxonomy_identifier, version):
        if (interactor_a in network and interactor_b not in network):
            nodes_to_add.add(interactor_b)

        elif (interactor_a not in network and interactor_b in network):
            nodes_to_add.add(interactor_a)

    network.add_nodes_from(nodes_to_add)


def add_protein_protein_interactions_from_biogrid(
        network: nx.Graph,
        experimental_system: Optional[Container[str]] = None,
        experimental_system_type: Optional[Container[str]] = None,
        interaction_throughput: Optional[Container[str]] = None,
        multi_validated_physical: bool = False,
        taxonomy_identifier: int = 9606,
        version: Optional[tuple[int, int, int]] = None) -> None:
    """
    Adds protein-protein interactions from BioGRID to a protein-protein
    interaction network.

    Args:
        network: The protein-protein interaction network.
        experimental_system: The accepted experimental evidence codes. If none
            are specified, any is accepted.
        experimental_system_type: The accepted categories of experimental
            evidence. If none are specified, any is accepted.
        interaction_throughput:  The accepted levels of interaction throughput.
            If none are specified, any is accepted.
        multi-validated physical: If True, add only multi-validated physical
            interactions.
        taxonomy_identifier: The taxonomy identifier.
        version: The version of the BioGRID database, if not passed, the latest.
    """
    for interactor_a, interactor_b in biogrid.get_protein_protein_interactions(
            experimental_system, experimental_system_type,
            interaction_throughput, multi_validated_physical,
            taxonomy_identifier, version):
        if (interactor_a in network and interactor_b in network and
                interactor_a != interactor_b):
            network.add_edge(interactor_a, interactor_b)
            network.edges[interactor_a, interactor_b]["BioGRID"] = 1.0


def add_proteins_from_intact(network: nx.Graph,
                             interaction_detection_methods: Optional[
                                 Container[str]] = None,
                             interaction_types: Optional[Container[str]] = None,
                             mi_score: float = 0.0,
                             taxonomy_identifier: int = 9606) -> None:
    """
    Add proteins interacting with proteins in a protein-protein interaction
    network from IntAct to the network.

    Args:
        network: The protein-protein interaction network.
        interaction_detection_methods: The accepted PSI-MI terms for interaction
            detection method. If none are specified, any is accepted.
        interaction_types: The accepted PSI-MI terms for interaction type. If
            none are specified, any is accepted.
        mi_score: The PSI-MI score threshold.
        taxonomy_identifier: The taxonomy identifier.
    """
    nodes_to_add = set()
    for interactor_a, interactor_b, _ in intact.get_protein_protein_interactions(
            interaction_detection_methods, interaction_types, mi_score,
            taxonomy_identifier):
        if (interactor_a in network and interactor_b not in network):
            nodes_to_add.add(interactor_b)

        elif (interactor_a not in network and interactor_b in network):
            nodes_to_add.add(interactor_a)

    network.add_nodes_from(nodes_to_add)


def add_protein_protein_interactions_from_intact(
        network: nx.Graph,
        interaction_detection_methods: Optional[Container[str]] = None,
        interaction_types: Optional[Container[str]] = None,
        mi_score: float = 0.0,
        taxonomy_identifier: int = 9606) -> None:
    """
    Adds protein-protein interactions from IntAct to a protein-protein
    interaction network.

    Args:
        network: The protein-protein interaction network.
        interaction_detection_methods: The accepted PSI-MI terms for interaction
            detection method. If none are specified, any is accepted.
        interaction_types: The accepted PSI-MI terms for interaction type.
            If none are specified, any is accepted.
        mi_score: The PSI-MI score threshold.
        taxonomy_identifier: The taxonomy identifier.
    """
    for interactor_a, interactor_b, score in intact.get_protein_protein_interactions(
            interaction_detection_methods, interaction_types, mi_score,
            taxonomy_identifier):
        if (interactor_a in network and interactor_b in network and
                interactor_a != interactor_b):
            if network.has_edge(interactor_a, interactor_b):
                network.edges[interactor_a, interactor_b]["IntAct"] = max(
                    score, network.edges[interactor_a,
                                         interactor_b].get("IntAct", 0.0))
            else:
                network.add_edge(interactor_a, interactor_b)
                network.edges[interactor_a, interactor_b]["IntAct"] = score


def add_proteins_from_mint(network: nx.Graph,
                           interaction_detection_methods: Optional[
                               Container[str]] = None,
                           interaction_types: Optional[Container[str]] = None,
                           mi_score: float = 0.0,
                           taxonomy_identifier: int = 9606) -> None:
    """
    Add proteins interacting with proteins in a protein-protein interaction
    network from MINT to the network.

    Args:
        network: The protein-protein interaction network.
        interaction_detection_methods: The accepted PSI-MI terms for interaction
            detection method. If none are specified, any is accepted.
        interaction_types: The accepted PSI-MI terms for interaction type.
            If none are specified, any is accepted.
        mi_score: The PSI-MI score threshold.
        taxonomy_identifier: The taxonomy identifier.
    """
    nodes_to_add = set()
    for interactor_a, interactor_b, _ in mint.get_protein_protein_interactions(
            interaction_detection_methods, interaction_types, mi_score,
            taxonomy_identifier):
        if (interactor_a in network and interactor_b not in network):
            nodes_to_add.add(interactor_b)

        elif (interactor_a not in network and interactor_b in network):
            nodes_to_add.add(interactor_a)

    network.add_nodes_from(nodes_to_add)


def add_protein_protein_interactions_from_mint(
        network: nx.Graph,
        interaction_detection_methods: Optional[Container[str]] = None,
        interaction_types: Optional[Container[str]] = None,
        mi_score: float = 0.0,
        taxonomy_identifier: int = 9606) -> None:
    """
    Adds protein-protein interactions from MINT to a protein-protein interaction
    network.

    Args:
        network: The protein-protein interaction network.
        interaction_detection_methods: The accepted PSI-MI terms for interaction
            detection method. If none are specified, any is accepted.
        interaction_types: The accepted PSI-MI terms for interaction type.
            If none are specified, any is accepted.
        mi_score: The PSI-MI score threshold.
        taxonomy_identifier: The taxonomy identifier.
    """
    for interactor_a, interactor_b, score in mint.get_protein_protein_interactions(
            interaction_detection_methods, interaction_types, mi_score,
            taxonomy_identifier):
        if (interactor_a in network and interactor_b in network and
                interactor_a != interactor_b):
            if network.has_edge(interactor_a, interactor_b):
                network.edges[interactor_a, interactor_b]["MINT"] = max(
                    score, network.edges[interactor_a,
                                         interactor_b].get("MINT", 0.0))
            else:
                network.add_edge(interactor_a, interactor_b)
                network.edges[interactor_a, interactor_b]["MINT"] = score


def add_proteins_from_reactome(
        network: nx.Graph,
        interaction_type: Optional[Container[str]] = None,
        interaction_context: Optional[Container[str]] = None,
        taxonomy_identifier: int = 9606) -> None:
    """
    Add proteins interacting with proteins in a protein-protein interaction
    network from Reactome to the network.

    Args:
        network: The protein-protein interaction network.
        interaction_type: The accepted interaction type annotation.
            If none are specified, any is accepted.
        interaction_context: The accepted interaction context annotation.
            If none are specified, any is accepted.
        taxonomy_identifier: The taxonomy identifier.
    """
    nodes_to_add = set()

    for interactor_a, interactor_b in reactome.get_protein_protein_interactions(
            interaction_type, interaction_context, taxonomy_identifier):
        if (interactor_a in network and interactor_b not in network):
            nodes_to_add.add(interactor_b)

        elif (interactor_a not in network and interactor_b in network):
            nodes_to_add.add(interactor_a)

    network.add_nodes_from(nodes_to_add)


def add_protein_protein_interactions_from_reactome(
    network: nx.Graph,
    interaction_type: Optional[Container[str]] = None,
    interaction_context: Optional[Container[str]] = None,
    taxonomy_identifier: int = 9606,
) -> None:
    """
    Adds protein-protein interactions from Reactome to a protein-protein
    interaction network.

    Args:
        network: The protein-protein interaction network.
        interaction_type: The accepted interaction type annotation.
            If none are specified, any is accepted.
        interaction_context: The accepted interaction context annotation.
            If none are specified, any is accepted.
        taxonomy_identifier: The taxonomy identifier.
    """
    for interactor_a, interactor_b in reactome.get_protein_protein_interactions(
            interaction_type, interaction_context, taxonomy_identifier):
        if (interactor_a in network and interactor_b in network and
                interactor_a != interactor_b):
            network.add_edge(interactor_a, interactor_b)
            network.edges[interactor_a, interactor_b]["Reactome"] = 1.0


def add_proteins_from_string(network: nx.Graph,
                             neighborhood: float = 0.0,
                             neighborhood_transferred: float = 0.0,
                             fusion: float = 0.0,
                             cooccurence: float = 0.0,
                             homology: float = 0.0,
                             coexpression: float = 0.0,
                             coexpression_transferred: float = 0.0,
                             experiments: float = 0.0,
                             experiments_transferred: float = 0.0,
                             database: float = 0.0,
                             database_transferred: float = 0.0,
                             textmining: float = 0.0,
                             textmining_transferred: float = 0.0,
                             combined_score: float = 0.0,
                             physical: bool = False,
                             taxonomy_identifier: int = 9606,
                             version: float = 11.5):
    """
    Add proteins interacting with proteins in a protein-protein interaction
    network from MINT to the network.

    Args:
        network: The protein-protein interaction network.
        neighborhood: The normal gene neighborhood score threshold.
        neighborhood_transferred: The transferred gene neighborhood score
            threshold.
        fusion: The gene fusion score threshold.
        cooccurrence: The gene cooccurrence score threshold.
        homology: The homology score threshold.
        coexpression: The normal coexpression score threshold.
        coexpression_transferred: The transferred coexpression score threshold.
        experiments: The normal experiments score threshold.
        experiments_transferred: The transferred experiments score threshold.
        database: The normal database score threshold.
        database_transferred: The transferred database score threshold.
        textmining: The normal textmining score threshold.
        textmining_transferred: The transferred textmining score threshold.
        combined_score: The combined score threshold.
        physical: If True, consider only physical interactions.
        taxonomy_identifier: The taxonomy identifier.
        version: The version of the STRING database.
    """
    nodes_to_add = set()

    for interactor_a, interactor_b, _ in string.get_protein_protein_interactions(
            neighborhood, neighborhood_transferred, fusion, cooccurence,
            homology, coexpression, coexpression_transferred, experiments,
            experiments_transferred, database, database_transferred, textmining,
            textmining_transferred, combined_score, physical,
            taxonomy_identifier, version):
        if (interactor_a in network and interactor_b not in network):
            nodes_to_add.add(interactor_b)

        elif (interactor_a not in network and interactor_b in network):
            nodes_to_add.add(interactor_a)

    network.add_nodes_from(nodes_to_add)


def add_protein_protein_interactions_from_string(
        network: nx.Graph,
        neighborhood: float = 0.0,
        neighborhood_transferred: float = 0.0,
        fusion: float = 0.0,
        cooccurence: float = 0.0,
        homology: float = 0.0,
        coexpression: float = 0.0,
        coexpression_transferred: float = 0.0,
        experiments: float = 0.0,
        experiments_transferred: float = 0.0,
        database: float = 0.0,
        database_transferred: float = 0.0,
        textmining: float = 0.0,
        textmining_transferred: float = 0.0,
        combined_score: float = 0.0,
        physical: bool = False,
        taxonomy_identifier: int = 9606,
        version: float = 11.5):
    """
    Adds protein-protein interactions from STRING to a protein-protein
    interaction network.

    Args:
        network: The protein-protein interaction network.
        neighborhood: The normal gene neighborhood score threshold.
        neighborhood_transferred: The transferred gene neighborhood score
            threshold.
        fusion: The gene fusion score threshold.
        cooccurrence: The gene cooccurrence score threshold.
        homology: The homology score threshold.
        coexpression: The normal coexpression score threshold.
        coexpression_transferred: The transferred coexpression score threshold.
        experiments: The normal experiments score threshold.
        experiments_transferred: The transferred experiments score threshold.
        database: The normal database score threshold.
        database_transferred: The transferred database score threshold.
        textmining: The normal textmining score threshold.
        textmining_transferred: The transferred textmining score threshold.
        combined_score: The combined score threshold.
        physical: If True, add only physical interactions.
        taxonomy_identifier: The taxonomy identifier.
        version: The version of the STRING database.
    """
    for interactor_a, interactor_b, score in string.get_protein_protein_interactions(
            neighborhood, neighborhood_transferred, fusion, cooccurence,
            homology, coexpression, coexpression_transferred, experiments,
            experiments_transferred, database, database_transferred, textmining,
            textmining_transferred, combined_score, physical,
            taxonomy_identifier, version):
        if (interactor_a in network and interactor_b in network and
                interactor_a != interactor_b):
            if network.has_edge(interactor_a, interactor_b):
                network.edges[interactor_a, interactor_b]["STRING"] = max(
                    score,
                    network.edges[interactor_a,
                                  interactor_b].get("STRING", 0.0),
                )
            else:
                network.add_edge(interactor_a, interactor_b)
                network.edges[interactor_a, interactor_b]["STRING"] = score


def get_databases(network: nx.Graph) -> tuple[str, ...]:
    """
    Returns the protein-protein interaction databases represented in a
    protein-protein interaction network.

    Args:
        network: The protein-protein interaction network.

    Returns:
        The databases of edge scores associated with edges in the
        protein-protein interaction network.
    """
    return tuple(
        sorted(
            set(database for edge in network.edges()
                for database in network.edges[edge]).intersection(
                    {"BioGRID", "IntAct", "MINT", "Reactome", "STRING"})))


def set_edge_weights(
    network: nx.Graph,
    weight: Callable[[dict[str, float]], float] = lambda confidence_scores: int(
        bool(confidence_scores.values())),
    attribute: str = "weight",
) -> None:
    """
    Set combined edge weights of a protein-protein interaction network.

    Args:
        network: The protein-protein interaction network.
        weight: The function to derive a combined confidence score from
            database-specific scores.
        attribute: The attribute name of the weight.
    """
    databases = get_databases(network)

    for edge in network.edges:
        network.edges[edge][attribute] = weight({
            database: network.edges[edge][database]
            for database in network.edges[edge]
            if database in databases
        })


def remove_edge_weights(network: nx.Graph, attribute: str = "weight") -> None:
    """
    Remove edge weights of a protein-protein interaction network.

    Args:
        network: The protein-protein interaction network.
        attribute: The attribute name of the weight.
    """
    for _, _, data in network.edges(data=True):
        if attribute in data:
            del data[attribute]


def get_modules(network: nx.Graph,
                module_size: int,
                module_size_combination: Callable[[list[int]],
                                                  float] = statistics.mean,
                algorithm: Callable[
                    [nx.Graph], list[set[Hashable]]] = modularization.louvain,
                resolution: float = 1.0,
                weight: str = "weight") -> list[nx.Graph]:
    """
    Returns modules of a protein-protein interaction network.

    Args:
        network: The protein-protein interaction network.
        module_size: The maximum module size. Community detection continues.
            subdivision of modles iteratively until it is met.
        module_size_combination: The function to derive a representative value
            from module sizes.
        algorithm: The community detection algorithm.
        resolution: The resolution parameter for modularity.
        weight: The attribute name of the edge weight.

    Returns:
        Modules of the protein-protein interaction network.
    """
    copied_network = network.copy()
    copied_network.remove_nodes_from(list(nx.isolates(copied_network)))

    communities = algorithm(copied_network, resolution, weight)

    while (module_size_combination(len(community) for community in communities)
           > module_size):

        subdivision = False
        for i, subdivided_community in enumerate(
                algorithm(copied_network.subgraph(communities[j]), resolution,
                          weight) for j in range(len(communities))):
            if len(subdivided_community) > 1:
                subdivision = True
                communities[i:i + 1] = subdivided_community

        if not subdivision:
            break

    return [network.subgraph(community) for community in communities]


def get_changes(
    network: nx.Graph,
    time: int,
    modification: str,
    site_combination: Optional[Callable[[Collection[float]], float]] = None
) -> tuple[float, ...]:
    """
    Returns the change distribution for a particular type of post-translational
    modification at a particular time of measurement.

    Args:
        network: The protein-protein interaction network.
        time: The time of measurement.
        modification: The type of post-translational modification.
        site_combination: An optional function to derive protein-specific
            changes from site-specific changes.

    Returns:
        Changes for a particular type of post-translational
        modification at a particular time of measurement.
    """
    changes = []
    for protein in network:
        sites = [
            network.nodes[protein][change]
            for change in network.nodes[protein]
            if len(change.split(" ")) == 3 and change.split(" ")[0] == str(time)
            and change.split(" ")[1] == modification
        ]

        if sites:
            if site_combination:
                changes.append(site_combination(sites))
            else:
                changes.extend(sites)

    return tuple(changes)


def convert_change_to_standard_score(
    network: nx.Graph,
    change: float,
    time: int,
    modification: str,
    site_combination: Optional[Callable[[Collection[float]], float]] = None
) -> float:
    """
    Converts a change to its corresponding modification- and time-specific
    standard score.

    Args:
        network: The protein-protein interaction network.
        change: The change to convert.
        time: The time of measurement.
        modification: The type of post-translational modification.
        site_combination: An optional function to derive protein-specific
            changes from site-specific changes.

    Returns:
        The standard score corresponding to the change.
    """
    changes = get_changes(network, time, modification, site_combination)
    mean = statistics.mean(changes)
    stdev = statistics.stdev(changes, xbar=mean)
    return (change - mean) / stdev


def convert_standard_score_to_change(
    network: nx.Graph,
    standard_score: float,
    time: int,
    modification: str,
    site_combination: Optional[Callable[[Collection[float]], float]] = None
) -> float:
    """
    Converts a modification- and time-specific standard score to its
    corresponding change.

    Args:
        network: The protein-protein interaction network.
        standard_score: The standard score to convert.
        time: The time of measurement.
        modification: The type of post-translational modification.
        site_combination: An optional function to derive protein-specific
            changes from site-specific changes.

    Returns:
        The change corresponding to the standard score.
    """
    changes = get_changes(network, time, modification, site_combination)
    mean = statistics.mean(changes)
    stdev = statistics.stdev(changes, xbar=mean)
    return standard_score * stdev + mean


def convert_change_to_quantile(
    network: nx.Graph,
    change: float,
    time: int,
    modification: str,
    site_combination: Optional[Callable[[Collection[float]], float]] = None
) -> float:
    """
    Converts a change to its corresponding modification- and time-specific
    quantile.

    Args:
        network: The protein-protein interaction network.
        change: The change to convert.
        time: The time of measurement.
        modification: The type of post-translational modification.
        site_combination: An optional function to derive protein-specific
            changes from site-specific changes.

    Returns:
        The quantile corresponding to the change.
    """
    changes = get_changes(network, time, modification, site_combination)
    return len([c for c in changes if c <= change]) / len(changes)


def convert_quantile_to_change(
    network: nx.Graph,
    quantile: float,
    time: int,
    modification: str,
    site_combination: Optional[Callable[[Collection[float]], float]] = None
) -> float:
    """
    Convert a modification- and time-specific quantile to its corresponding
    change.

    Args:
        network: The protein-protein interaction network.
        quantile: The quantile to convert.
        time: The time of measurement.
        modification: The type of post-translational modification.
        site_combination: An optional function to derive protein-specific
            changes from site-specific changes.

    Returns:
        The change corresponding to the quantile.
    """
    changes = get_changes(network, time, modification, site_combination)
    return sorted(changes)[math.floor(quantile * (len(changes) - 1))]


def get_change_enrichment(
    network: nx.Graph,
    modules: list[nx.Graph],
    changes: tuple[float, float] = (-1.0, 1.0),
    convert_change: Callable[
        [nx.Graph, float, int, str, Callable[[Collection[float]],
                                             float]], float] = lambda network,
    change, time, modification, site_combination: change,
    site_combination: Optional[Callable[[Collection[float]], float]] = None,
    enrichment_test: Callable[[int, int, int, int],
                              float] = test.hypergeometric,
    multiple_testing_correction: Callable[
        [dict[tuple[nx.Graph, int, str], float]],
        dict[tuple[nx.Graph, int, str], float]] = correction.benjamini_hochberg,
) -> dict[nx.Graph, dict[int, dict[str, float]]]:
    """
    Test modules for enrichment of large protein-specific changes for each time
    of measurement and type of post-translational modification with respect to
    the protein-protein interaction network.

    Args:
        network: The protein-protein interaction network.
        modules: The modules of the protein-protein interaction network.
        changes: Proteins are classified by whether their representative change
            exceeds this range.
        convert_change: The function used convert the bounds to log2-fold
            changes.
        site_combination: An optional function to derive protein-specific
            changes from site-specific changes.
        enrichment_test: The statistical test used to assess enrichment of
            proteins with large changes by a module.
        multiple_testing_correction: The procedure to correct for multiple
            testing of multiple modules, times of measurement and types of
            post-translational modification.

    Returns:
        Adjusted p-values for the enrichment of large protein-specific changes
        by each module for each time of measurement and type of
        post-translational modification.
    """
    p_values = {}
    for time in get_times(network):
        for modification in get_post_translational_modifications(network, time):
            change_range = (convert_change(network, changes[0], time,
                                           modification, site_combination),
                            convert_change(network, changes[1], time,
                                           modification, site_combination))

            modified_proteins = len([
                change for change in get_changes(network, time, modification,
                                                 site_combination) if change
            ])

            modified_module_proteins = [
                len([
                    change for change in get_changes(module, time, modification,
                                                     site_combination) if change
                ]) for module in modules
            ]

            target_proteins = len([
                change for change in get_changes(network, time, modification,
                                                 site_combination)
                if change <= change_range[0] or change >= change_range[1]
            ])

            target_module_proteins = [
                len([
                    change
                    for change in get_changes(module, time, modification,
                                              site_combination)
                    if change <= change_range[0] or change >= change_range[1]
                ])
                for module in modules
            ]

            p_values.update({(module, time, modification):
                             enrichment_test(target_module_proteins[i],
                                             modified_proteins, target_proteins,
                                             modified_module_proteins[i])
                             for i, module in enumerate(modules)})

    p_values = multiple_testing_correction(p_values)

    return {
        module: {
            time: {
                modification: p_values[(module, time, modification)]
                for modification in get_post_translational_modifications(
                    network, time)
            } for time in get_times(network)
        } for module in modules
    }


def get_change_location(
    network: nx.Graph,
    modules: list[nx.Graph],
    site_combination: Optional[Callable[[Collection[float]], float]] = None,
    location_test: Callable[[Collection[float], Collection[float]],
                            float] = test.wilcoxon,
    multiple_testing_correction: Callable[
        [dict[tuple[nx.Graph, int, str], float]],
        dict[tuple[nx.Graph, int, str], float]] = correction.benjamini_hochberg,
) -> dict[nx.Graph, dict[int, dict[str, float]]]:
    """
    Test modules for difference tendencies in protein-specific changes for each
    time of measurement and type of post-translational modification with respect
    to the remaining protein-protein interaction network.

    Args:
        network: The protein-protein interaction network.
        modules: The modules of the protein-protein interaction network.
        site_combination: An optional function to derive protein-specific
            changes from site-specific changes.
        location_test: The statistical test used to assess location of
            proteins with large changes by a module.
        multiple_testing_correction: The procedure to correct for multiple
            testing of multiple modules, times of measurement and types of
            post-translational modification.

    Returns:
        Adjusted p-values for the difference in central tendencies in
        protein-specific changes of modules for each time of measurement and
        type of post-translational modification.
    """
    p_values = {}
    for time in get_times(network):
        for modification in get_post_translational_modifications(network, time):
            network_changes = [
                get_changes(nx.union_all([m
                                          for m in modules
                                          if m != module]), time, modification,
                            site_combination)
                for module in modules
            ]

            module_changes = [
                get_changes(module, time, modification, site_combination)
                for module in modules
            ]

            p_values.update({(module, time, modification):
                             location_test(module_changes[i],
                                           network_changes[i])
                             for i, module in enumerate(modules)})

    p_values = multiple_testing_correction(p_values)

    return {
        module: {
            time: {
                modification: p_values[(module, time, modification)]
                for modification in get_post_translational_modifications(
                    network, time)
            } for time in get_times(network)
        } for module in modules
    }


def get_gene_ontology_enrichment(
    networks: list[nx.Graph],
    enrichment_test: Callable[[int, int, int, int],
                              float] = test.hypergeometric,
    multiple_testing_correction: Callable[[dict[tuple[
        nx.Graph, str], float]], dict[tuple[nx.Graph, str],
                                      float]] = correction.benjamini_hochberg,
    taxonomy_identifier: int = 9606,
    namespaces: Container[str] = ("cellular_component", "molecular_function",
                                  "biological_process"),
    annotation_as_reference: bool = True
) -> dict[nx.Graph, dict[tuple[str, str], float]]:
    """
    Test the protein-protein interaction networks for enrichment of Gene
    Ontology terms.

    Args:
        networks: The protein-protein interaction networks.
        enrichment_test: The statistical test used to assess enrichment of a
            Gene Ontology term.
        multiple_testing_correction: The procedure to correct for multiple
            testing of multiple terms and networks.
        taxonomy_identifier: The taxonomy identifier.
        namespaces: The Gene Ontology namespaces.
        annotation_as_reference: If True, compute enrichment with respect to the
            entire species-specific Gene Ontology annotation in namespaces, else
            with respect to the union of the protein-protein interaction
            networks.

    Returns:
        Adjusted p-values for the enrichment of each Gene Ontology term by each
        network.
    """
    name, go_id = {}, {}
    for term in gene_ontology.get_ontology(namespaces):
        name[term["id"]] = term["name"]
        for alt_id in term["alt_id"]:
            if alt_id not in go_id:
                go_id[alt_id] = set()
            go_id[alt_id].add(term["id"])

    annotation = {}
    for protein, term in gene_ontology.get_annotation(
            taxonomy_identifier, gene_ontology.convert_namespaces(namespaces)):
        if annotation_as_reference or any(
                protein in network.nodes() for network in networks):
            for primary_term in go_id.get(term, {term}):
                if primary_term not in annotation:
                    annotation[primary_term] = set()
                annotation[primary_term].add(protein)

    annotation = {
        term: proteins for term, proteins in annotation.items() if proteins
    }

    annotated_proteins = set.union(*annotation.values())

    annotated_network_proteins = {
        network: len(annotated_proteins.intersection(network.nodes()))
        for network in networks
    }

    intersection = {
        network: {
            term: len(annotation[term].intersection(network.nodes()))
            for term in annotation
        } for network in networks
    }

    p_value = multiple_testing_correction({
        (network, term): enrichment_test(intersection[network][term],
                                         len(annotated_proteins),
                                         len(annotation[term]),
                                         annotated_network_proteins[network])
        for term in annotation for network in networks
    })

    return {
        network: {(term, name[term]): p_value[(network, term)]
                  for term in annotation} for network in networks
    }


def export(network: nx.Graph, basename: str) -> None:
    """
    Exports the protein-protein interaction network.

    Args:
        network: The protein-protein interaction network.
        basename: The base file name.
    """
    nx.write_graphml_xml(network,
                         "{}.graphml".format(basename),
                         named_key_ids=True,
                         infer_numeric_types=True)
