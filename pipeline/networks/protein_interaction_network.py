"""protein-protein interaction network"""
import bisect
import math
import os
import re
import statistics
from typing import (Callable, Collection, Container, Hashable, Iterable,
                    MutableSequence, Optional)

import networkx as nx
import pandas as pd
import scipy.stats
from analysis import correction, modularization
from databases import biogrid, corum, intact, mint, reactome, string, uniprot


def get_network() -> nx.Graph:
    """
    Returns an initialized protein-protein interaction network.

    Returns
        The protein-protein interaction network.

    """
    return nx.Graph()


def add_genes_from_table(
    network: nx.Graph,
    file_name: str,
    gene_accession_column: int | str,
    gene_accession_format: re.Pattern = re.compile("^(.+?)$"),
    sheet_name: int | str = 0,
    header: int = 0,
) -> None:
    """
    Parse UniProt gene accessions from a tabular file and add the corresponding
    primary UniProt protein accessions to a protein-protein interaction network.

    Args:
        network: The protein-protein interaction network.
        file_name: The file location of the file.
        gene_accession_column: The column containing UniProt gene accessions.
        gene_accession_format: A regular expression to extract gene accessions
            from an entry.
        sheet_name: The sheet to parse gene accessions from.
        header: The index of the header row.
        organism: The NCBI taxonomy identifier for the organism of interest.
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

    network.add_nodes_from(genes)


def map_genes(network: nx.Graph, organism: int = 9606) -> None:
    """
    Map genes in a protein-protein interaction network to their primary
    UniProt protein identifiers and remove genes not present in Swiss-Prot.

    Args:
        network: The protein-protein interaction network to map genes from.
        organism: The NCBI taxonomy identifier for the organism of interest.
    """
    mapping, gene_name, protein_name = {}, {}, {}
    for accessions, gene, protein in uniprot.get_swiss_prot_entries(organism):
        if gene in network.nodes():
            mapping[gene] = accessions[0]
            gene_name[accessions[0]] = gene
            protein_name[accessions[0]] = protein

    network.remove_nodes_from(set(network.nodes()).difference(mapping))
    nx.relabel_nodes(network, mapping, copy=False)

    for protein in network.nodes():
        network.nodes[protein]["gene"] = gene_name[protein]
        network.nodes[protein]["protein"] = protein_name[protein]


def add_proteins_from_table(
    network: nx.Graph,
    file_name: str,
    protein_accession_column: int | str,
    protein_accession_format: re.Pattern = re.compile("^(.+?)$"),
    time: int = 0,
    modification: str = "",
    position_column: int | str = "",
    position_format: re.Pattern = re.compile("^(.+?)$"),
    replicates: Optional[Collection[int | str]] = None,
    sheet_name: int | str = 0,
    header: int = 0,
    number_sites: int = 100,
    number_replicates: int = 1,
    replicate_combination: Callable[[Iterable[float]], float] = statistics.mean,
    measurement_conversion: Callable[
        [Iterable[float], Callable[[Iterable[float]], float]],
        float] = lambda measurements, combination: math.log2(
            combination(measurements))
) -> None:
    """
    Parse UniProt protein accessions and measurements of measurements from a
    tabular file and add the corresponding primary UniProt protein accessions to
    a protein-protein interaction network.

    Args:
        network: The protein-protein interaction network.
        file_name: The file location of the file.
        protein_accession_column: The column containing UniProt protein
            accessions.
        protein_accession_format: A regular expression to extract protein
            accessions from a corresponding entry.
        time: The time of measurement to associate with measurements.
        modification: An identifier for the type of post-translational
            modification to associate with measurements.
        position_column: The column containing sites corresponding to
            measurements.
        replicates: The columns containing replicates of measurements.
        sheet_name: The sheet to parse protein accessions from.
        header: The index of the header row.
        number_sites: The maximum number of measurements to associate with a
            protein-accession, prioritized by largest absolute value.
        number_replicates: The minimum number of replicates to accept a
            measurement.
        replicate_combination: A function to combine replicates into a single
            measurement.
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
            ([position_column] if position_column else []) + list(replicates),
            dtype={
                protein_accession_column: str,
                **{
                    position_column: str if position_column else {}
                },
                **{column: float for column in replicates},
            },
        )
    else:
        table = pd.read_table(
            file_name,
            header=header,
            usecols=[protein_accession_column] +
            ([position_column] if position_column else []) + list(replicates),
            dtype={
                protein_accession_column: str,
                **{
                    position_column: str if position_column else {}
                },
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

            if len(measurements) >= min(number_replicates, len(replicates)):
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

    network.add_nodes_from(proteins)
    for protein, sites in proteins.items():
        measurements = sorted(
            sorted(
                sites,
                key=lambda item: abs(item[1]),
            )[-number_sites:])

        for i, (_, measurement) in enumerate(measurements):
            network.nodes[protein][
                f"{time} {modification} {i + 1}"] = measurement


def map_proteins(network: nx.Graph, organism: int = 9606) -> None:
    """
    Map proteins in a protein-protein interaction network to their primary
    UniProt identifiers and remove proteins not present in Swiss-Prot. Isoform
    identifiers are maintained, but not transferred.

    Args:
        network: The protein-protein interaction network to map proteins from.
        organism: The NCBI taxonomy identifier for the organism of interest.
    """
    mapping, gene_name, protein_name = {}, {}, {}
    for accessions, gene, protein in uniprot.get_swiss_prot_entries(organism):
        for node in network.nodes():
            if node.split("-")[0] in accessions:
                if "-" in node and node.split("-")[1].isnumeric(
                ) and accessions.index(node.split("-")[0]) == 0:
                    mapping[node] = node
                    gene_name[node] = gene
                    protein_name[node] = protein
                else:
                    mapping[node] = accessions[0]
                    gene_name[accessions[0]] = gene
                    protein_name[accessions[0]] = protein

    network.remove_nodes_from(set(network.nodes()).difference(mapping))
    nx.relabel_nodes(network, mapping, copy=False)

    for protein in network.nodes():
        network.nodes[protein]["gene"] = gene_name[protein]
        network.nodes[protein]["protein"] = protein_name[protein]


def get_proteins(
    network: nx.Graph,
    time: int,
    modification: str,
    site_combination: Callable[[Collection[float]],
                               float] = lambda sites: max(sites, key=abs),
    combined_measurement_filter: Callable[[float],
                                          bool] = bool) -> frozenset[str]:
    """
    Returns proteins of a protein-protein interaction network with specified
    measurements.

    Args:
        network: The protein-protein interaction network.
        time: The time of measurement.
        modification: The type of post-translational modification.
        site_combination: An optional function to derive protein-specific
            measurements from site-specific measurements.
        combined_measurement_filter: The predicate to filter combined
            measurements.

    Returns:
        Proteins whose combined measurement of a particular type of
        post-translational modification at a particular time of measurement
        fulfil the predicate.
    """
    proteins = []
    for protein in network:
        sites = [
            network.nodes[protein][measurement]
            for measurement in network.nodes[protein]
            if len(measurement.split(" ")) == 3 and measurement.split(" ")[0] ==
            str(time) and measurement.split(" ")[1] == modification
        ]

        if sites and combined_measurement_filter(site_combination(sites)):
            proteins.append(protein)

    return frozenset(proteins)


def get_times(network: nx.Graph) -> tuple[int, ...]:
    """
    Returns the times of measurement represented in a protein-protein
    interaction network.

    Args:
        network: The protein-protein interaction network.

    Returns:
        The times of measurement associated with any measurements of any protein
        in the protein-protein interaction network.
    """
    return tuple(
        sorted(
            set(
                int(measurement.split(" ")[0])
                for protein in network
                for measurement in network.nodes[protein]
                if len(measurement.split(" ")) == 3 and
                measurement.split(" ")[0].isnumeric())))


def get_modifications(network: nx.Graph, time: int) -> tuple[str, ...]:
    """
    Returns the types of post-translational modification represented in a
    protein-protein interaction network at a particular time of measurement.

    Args:
        network: The protein-protein interaction network.
        time: The time of measurement.

    Returns:
        The types of post-translational modification associated with any
        measurements of any protein in the protein-protein interaction network
        at a particular time of measurement.
    """
    return tuple(
        sorted(
            set(
                measurement.split(" ")[1]
                for protein in network
                for measurement in network.nodes[protein]
                if len(measurement.split(" ")) == 3 and
                measurement.split(" ")[0] == str(time))))


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
        int(measurement.split(" ")[2])
        for protein in network
        for measurement in network.nodes[protein]
        if len(measurement.split(" ")) == 3 and measurement.split(" ")[0] ==
        str(time) and measurement.split(" ")[1] == modification)


def set_post_translational_modification(network: nx.Graph) -> None:
    """
    Annotate proteins with summary of type of corresponding modifications.

    Args:
        network: The protein-protein interaction network.
    """
    for time in get_times(network):
        for protein in network:
            network.nodes[protein][
                f"post-translational modification {time}"] = " ".join(
                    sorted(
                        set(
                            measurement.split(" ")[1]
                            for measurement in network.nodes[protein]
                            if len(measurement.split(" ")) == 3 and
                            measurement.split(" ")[0] == str(time))))


def set_measurements(
    network: nx.Graph,
    site_combination: Callable[[Collection[float]],
                               float] = lambda sites: max(sites, key=abs),
    measurements: tuple[float, float] = (-1.0, 1.0),
    convert_measurement: Callable[
        [float, Collection[float]],
        float] = lambda measurement, measurements: measurement
) -> None:
    """
    Annotate nodes with summary of measurements.

    Args:
        network: The protein-protein interaction network.
        site_combination: An optional function to derive protein-specific
            measurements from site-specific measurements.
        measurements: Proteins are categorized by whether their representative
            exceed either this range, the range defined by half the bounds or
            none.
        convert_measurement: The function used convert the bounds to log2-fold
            measurements.
    """
    times = get_times(network)
    modifications = {time: get_modifications(network, time) for time in times}

    measurement_range = {
        time: {
            modification: (convert_measurement(
                measurements[0],
                get_measurements(network, time, modification, site_combination),
            ),
                           convert_measurement(
                               measurements[1],
                               get_measurements(network, time, modification,
                                                site_combination)))
            for modification in modifications[time]
        } for time in times
    }

    for time in times:
        for protein in network:
            classification = {}
            for modification in modifications[time]:
                sites = [
                    network.nodes[protein][measurement]
                    for measurement in network.nodes[protein]
                    if len(measurement.split(" ")) == 3 and
                    measurement.split(" ")[0] == str(time) and
                    measurement.split(" ")[1] == modification
                ]

                if sites:
                    combined_sites = site_combination(sites)

                    if combined_sites >= 1.0 * measurement_range[time][
                            modification][1]:
                        classification[modification] = "up"
                    elif 1.0 * measurement_range[time][modification][
                            1] > combined_sites >= 0.5 * measurement_range[
                                time][modification][1]:
                        classification[modification] = "mid up"
                    elif 0.5 * measurement_range[time][modification][
                            1] > combined_sites > 0.5 * measurement_range[time][
                                modification][0]:
                        classification[modification] = "mid"
                    elif 0.5 * measurement_range[time][modification][
                            0] >= combined_sites > 1.0 * measurement_range[
                                time][modification][0]:
                        classification[modification] = "mid down"
                    else:
                        classification[modification] = "down"

            if classification:
                if set(classification.values()) == {"up"}:
                    network.nodes[protein][f"measurement {time}"] = "up"
                elif set(classification.values()) == {"mid up", "up"}:
                    network.nodes[protein][f"measurement {time}"] = "up"
                elif set(classification.values()) == {"mid", "mid up", "up"}:
                    network.nodes[protein][f"measurement {time}"] = "up"

                elif set(classification.values()) == {"mid up"}:
                    network.nodes[protein][f"measurement {time}"] = "mid up"
                elif set(classification.values()) == {"mid", "mid up"}:
                    network.nodes[protein][f"measurement {time}"] = "mid up"

                elif set(classification.values()) == {"mid"}:
                    network.nodes[protein][f"measurement {time}"] = "mid"

                elif set(classification.values()) == {"mid", "mid down"}:
                    network.nodes[protein][f"measurement {time}"] = "mid down"
                elif set(classification.values()) == {"mid down"}:
                    network.nodes[protein][f"measurement {time}"] = "mid down"

                elif set(
                        classification.values()) == {"mid", "mid down", "down"}:
                    network.nodes[protein][f"measurement {time}"] = "down"
                elif set(classification.values()) == {"mid down", "down"}:
                    network.nodes[protein][f"measurement {time}"] = "down"
                elif set(classification.values()) == {"down"}:
                    network.nodes[protein][f"measurement {time}"] = "down"
                else:
                    network.nodes[protein][f"measurement {time}"] = " ".join(
                        f"{modification} {classification[modification]}"
                        for modification in sorted(classification))

            else:
                network.nodes[protein][f"measurement {time}"] = "mid"


def get_neighbors_from_biogrid(
        network: nx.Graph,
        experimental_system: Optional[Container[str]] = None,
        experimental_system_type: Optional[Container[str]] = None,
        interaction_throughput: Optional[Container[str]] = None,
        multi_validated_physical: bool = False,
        organism: int = 9606,
        version: Optional[tuple[int, int, int]] = None) -> set[str]:
    """
    Returns proteins interacting with proteins in a protein-protein interaction
    network from BioGRID.

    Args:
        network: The protein-protein interaction network.
        experimental_system: The accepted experimental evidence codes. If none
            are specified, any are accepted.
        experimental_system_type: The accepted categories of experimental
            evidence. If none are specified, any are accepted.
        interaction_throughput:  The accepted levels of interaction throughput.
            If none are specified, any are accepted.
        multi-validated physical: If True, consider only multi-validated
            physical interactions.
        organism: The NCBI taxonomy identifier for the organism of interest.
        version: The version of the BioGRID database, if not specified or not
            consisting of three entries, the latest.

    Returns:
        Neighbors of the protein-protein interaction network in BioGRID.
    """
    neighbors = set()
    for interactor_a, interactor_b in biogrid.get_protein_interactions(
            experimental_system, experimental_system_type,
            interaction_throughput, multi_validated_physical, organism,
            version):
        if (interactor_a in network and interactor_b not in network):
            neighbors.add(interactor_b)

        elif (interactor_a not in network and interactor_b in network):
            neighbors.add(interactor_a)

    return neighbors


def add_protein_interactions_from_biogrid(
        network: nx.Graph,
        experimental_system: Optional[Container[str]] = None,
        experimental_system_type: Optional[Container[str]] = None,
        interaction_throughput: Optional[Container[str]] = None,
        multi_validated_physical: bool = False,
        organism: int = 9606,
        version: Optional[tuple[int, int, int]] = None) -> None:
    """
    Adds protein-protein interactions from BioGRID to a protein-protein
    interaction network.

    Args:
        network: The protein-protein interaction network.
        experimental_system: The accepted experimental evidence codes. If none
            are specified, any are accepted.
        experimental_system_type: The accepted categories of experimental
            evidence. If none are specified, any are accepted.
        interaction_throughput:  The accepted levels of interaction throughput.
            If none are specified, any are accepted.
        multi-validated physical: If True, add only multi-validated physical
            interactions.
        organism: The NCBI taxonomy identifier for the organism of interest.
        version: The version of the BioGRID database, if not specified or not
            consisting of three entries, the latest.
    """
    for interactor_a, interactor_b in biogrid.get_protein_interactions(
            experimental_system, experimental_system_type,
            interaction_throughput, multi_validated_physical, organism,
            version):
        if (interactor_a in network and interactor_b in network and
                interactor_a != interactor_b):
            network.add_edge(interactor_a, interactor_b)
            network.edges[interactor_a, interactor_b]["BioGRID"] = 1.0


def get_neighbors_from_corum(network: nx.Graph,
                             purification_methods: Optional[
                                 Container[str]] = None,
                             organism: int = 9606) -> set[str]:
    """
    Returns proteins interacting with proteins in a protein-protein interaction
    network from CORUM.

    Args:
        network: The protein-protein interaction network.
        purification_methods: The accepted PSI-MI identifiers or terms for the
            protein complex purification method. If none are specified, any are
            accepted.
        organism: The NCBI taxonomy identifier for the organism of interest.

    Returns:
        Neighbors of the protein-protein interaction network in CORUM.
    """
    neighbors = set()
    for interactor_a, interactor_b in corum.get_protein_interactions(
            purification_methods, organism):
        if (interactor_a in network and interactor_b not in network):
            neighbors.add(interactor_b)

        elif (interactor_a not in network and interactor_b in network):
            neighbors.add(interactor_a)

    return neighbors


def add_protein_interactions_from_corum(network: nx.Graph,
                                        purification_methods: Optional[
                                            Container[str]] = None,
                                        organism: int = 9606) -> None:
    """
    Adds protein-protein interactions from CORUM to a protein-protein
    interaction network.

    Args:
        network: The protein-protein interaction network.
        purification_methods: The accepted PSI-MI identifiers or terms for the
            protein complex purification method. If none are specified, any are
            accepted.
        organism: The NCBI taxonomy identifier for the organism of interest.
    """
    for interactor_a, interactor_b in corum.get_protein_interactions(
            purification_methods, organism):
        if (interactor_a in network and interactor_b in network and
                interactor_a != interactor_b):
            network.add_edge(interactor_a, interactor_b)
            network.edges[interactor_a, interactor_b]["CORUM"] = 1.0


def get_neighbors_from_intact(network: nx.Graph,
                              interaction_detection_methods: Optional[
                                  Iterable[str]] = None,
                              interaction_types: Optional[Iterable[str]] = None,
                              psi_mi_score: float = 0.0,
                              organism: int = 9606) -> set[str]:
    """
    Returns proteins interacting with proteins in a protein-protein interaction
    network from IntAct to the network.

    Args:
        network: The protein-protein interaction network.
        interaction_detection_methods: The accepted PSI-MI identifiers or terms
            for interaction detection method. If none are specified, any are
            accepted.
        interaction_types: The accepted PSI-MI identifiers or terms for
            interaction type. If none are specified, any are accepted.
        psi_mi_score: The PSI-MI score threshold.
        organism: The NCBI taxonomy identifier for the organism of interest.

    Returns:
        Neighbors of the protein-protein interaction network in IntAct.
    """
    neighbors = set()
    for interactor_a, interactor_b, _ in intact.get_protein_interactions(
            interaction_detection_methods, interaction_types, psi_mi_score,
            organism):
        if (interactor_a in network and interactor_b not in network):
            neighbors.add(interactor_b)

        elif (interactor_a not in network and interactor_b in network):
            neighbors.add(interactor_a)

    return neighbors


def add_protein_interactions_from_intact(
        network: nx.Graph,
        interaction_detection_methods: Optional[Iterable[str]] = None,
        interaction_types: Optional[Iterable[str]] = None,
        psi_mi_score: float = 0.0,
        organism: int = 9606) -> None:
    """
    Adds protein-protein interactions from IntAct to a protein-protein
    interaction network.

    Args:
        network: The protein-protein interaction network.
        interaction_detection_methods: The accepted PSI-MI identifiers or terms
            for interaction detection method. If none are specified, any are
            accepted.
        interaction_types: The accepted PSI-MI identifiers or terms for
            interaction type. If none are specified, any are accepted.
        psi_mi_score: The PSI-MI score threshold.
        organism: The NCBI taxonomy identifier for the organism of interest.
    """
    for interactor_a, interactor_b, score in intact.get_protein_interactions(
            interaction_detection_methods, interaction_types, psi_mi_score,
            organism):
        if (interactor_a in network and interactor_b in network and
                interactor_a != interactor_b):
            if network.has_edge(interactor_a, interactor_b):
                network.edges[interactor_a, interactor_b]["IntAct"] = max(
                    score, network.edges[interactor_a,
                                         interactor_b].get("IntAct", 0.0))
            else:
                network.add_edge(interactor_a, interactor_b)
                network.edges[interactor_a, interactor_b]["IntAct"] = score


def get_neighbors_from_mint(network: nx.Graph,
                            interaction_detection_methods: Optional[
                                Iterable[str]] = None,
                            interaction_types: Optional[Iterable[str]] = None,
                            psi_mi_score: float = 0.0,
                            organism: int = 9606) -> set[str]:
    """
    Returns proteins interacting with proteins in a protein-protein interaction
    network from MINT.

    Args:
        network: The protein-protein interaction network.
        interaction_detection_methods: The accepted PSI-MI identifiers or terms
            for interaction detection method. If none are specified, any are
            accepted.
        interaction_types: The accepted PSI-MI identifiers or terms for
            interaction type. If none are specified, any are accepted.
        psi_mi_score: The PSI-MI score threshold.
        organism: The NCBI taxonomy identifier for the organism of interest.

    Returns:
        Neighbors of the protein-protein interaction network in MINT.
    """
    neighbors = set()
    for interactor_a, interactor_b, _ in mint.get_protein_interactions(
            interaction_detection_methods, interaction_types, psi_mi_score,
            organism):
        if (interactor_a in network and interactor_b not in network):
            neighbors.add(interactor_b)

        elif (interactor_a not in network and interactor_b in network):
            neighbors.add(interactor_a)

    return neighbors


def add_protein_interactions_from_mint(
        network: nx.Graph,
        interaction_detection_methods: Optional[Iterable[str]] = None,
        interaction_types: Optional[Iterable[str]] = None,
        psi_mi_score: float = 0.0,
        organism: int = 9606) -> None:
    """
    Adds protein-protein interactions from MINT to a protein-protein interaction
    network.

    Args:
        network: The protein-protein interaction network.
        interaction_detection_methods: The accepted PSI-MI identifiers or terms
            for interaction detection method. If none are specified, any are
            accepted.
        interaction_types: The accepted PSI-MI identifiers or terms for
            interaction type. If none are specified, any are accepted.
        psi_mi_score: The PSI-MI score threshold.
        organism: The NCBI taxonomy identifier for the organism of interest.
    """
    for interactor_a, interactor_b, score in mint.get_protein_interactions(
            interaction_detection_methods, interaction_types, psi_mi_score,
            organism):
        if (interactor_a in network and interactor_b in network and
                interactor_a != interactor_b):
            if network.has_edge(interactor_a, interactor_b):
                network.edges[interactor_a, interactor_b]["MINT"] = max(
                    score, network.edges[interactor_a,
                                         interactor_b].get("MINT", 0.0))
            else:
                network.add_edge(interactor_a, interactor_b)
                network.edges[interactor_a, interactor_b]["MINT"] = score


def get_neighbors_from_reactome(
        network: nx.Graph,
        interaction_type: Optional[Container[str]] = None,
        interaction_context: Optional[Container[str]] = None,
        organism: int = 9606) -> set[str]:
    """
    Returns proteins interacting with proteins in a protein-protein interaction
    network from Reactome.

    Args:
        network: The protein-protein interaction network.
        interaction_type: The accepted interaction type annotation.
            If none are specified, any are accepted.
        interaction_context: The accepted interaction context annotation.
            If none are specified, any are accepted.
        organism: The NCBI taxonomy identifier for the organism of interest.

    Returns:
        Neighbors of the protein-protein interaction network in Reactome.
    """
    neighbors = set()

    for interactor_a, interactor_b in reactome.get_protein_interactions(
            interaction_type, interaction_context, organism):
        if (interactor_a in network and interactor_b not in network):
            neighbors.add(interactor_b)

        elif (interactor_a not in network and interactor_b in network):
            neighbors.add(interactor_a)

    return neighbors


def add_protein_interactions_from_reactome(
    network: nx.Graph,
    interaction_type: Optional[Container[str]] = None,
    interaction_context: Optional[Container[str]] = None,
    organism: int = 9606,
) -> None:
    """
    Adds protein-protein interactions from Reactome to a protein-protein
    interaction network.

    Args:
        network: The protein-protein interaction network.
        interaction_type: The accepted interaction type annotation.
            If none are specified, any are accepted.
        interaction_context: The accepted interaction context annotation.
            If none are specified, any are accepted.
        organism: The NCBI taxonomy identifier for the organism of interest.
    """
    for interactor_a, interactor_b in reactome.get_protein_interactions(
            interaction_type, interaction_context, organism):
        if (interactor_a in network and interactor_b in network and
                interactor_a != interactor_b):
            network.add_edge(interactor_a, interactor_b)
            network.edges[interactor_a, interactor_b]["Reactome"] = 1.0


def get_neighbors_from_string(network: nx.Graph,
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
                              organism: int = 9606,
                              version: float = 11.5) -> set[str]:
    """
    Add proteins interacting with proteins in a protein-protein interaction
    network from STRING to the network.

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
        organism: The NCBI taxonomy identifier for the organism of interest.
        version: The version of the STRING database.

    Returns:
        Neighbors of the protein-protein interaction network in STRING.
    """
    neighbors = set()

    for interactor_a, interactor_b, _ in string.get_protein_interactions(
            neighborhood, neighborhood_transferred, fusion, cooccurence,
            homology, coexpression, coexpression_transferred, experiments,
            experiments_transferred, database, database_transferred, textmining,
            textmining_transferred, combined_score, physical, organism,
            version):
        if (interactor_a in network and interactor_b not in network):
            neighbors.add(interactor_b)

        elif (interactor_a not in network and interactor_b in network):
            neighbors.add(interactor_a)

    return neighbors


def add_protein_interactions_from_string(network: nx.Graph,
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
                                         organism: int = 9606,
                                         version: float = 11.5) -> None:
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
        organism: The NCBI taxonomy identifier for the organism of interest.
        version: The version of the STRING database.
    """
    for interactor_a, interactor_b, score in string.get_protein_interactions(
            neighborhood, neighborhood_transferred, fusion, cooccurence,
            homology, coexpression, coexpression_transferred, experiments,
            experiments_transferred, database, database_transferred, textmining,
            textmining_transferred, combined_score, physical, organism,
            version):
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
                    ("BioGRID", "CORUM", "IntAct", "MINT", "Reactome",
                     "STRING"))))


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


def get_communities(
        network: nx.Graph,
        community_size: int,
        community_size_combination: Callable[[Iterable[int]],
                                             float] = statistics.mean,
        algorithm: Callable[
            [nx.Graph, float, str],
            MutableSequence[set[Hashable]]] = modularization.louvain,
        resolution: float = 1.0,
        weight: str = "weight") -> tuple[nx.Graph, ...]:
    """
    Returns communities of a protein-protein interaction network.

    Args:
        network: The protein-protein interaction network.
        community_size: The maximum community size. Communities are iteratively
            subdivided until it is met.
        community_size_combination: The function to derive a decisive value to
            from meeting the threshhold from the community sizes.
        algorithm: The community detection algorithm.
        resolution: The resolution parameter for modularity.
        weight: The attribute name of the edge weight.

    Returns:
        Modules of the protein-protein interaction network.
    """
    if not network.number_of_edges():
        return tuple()

    copied_network = network.copy()
    copied_network.remove_nodes_from(list(nx.isolates(copied_network)))

    communities = algorithm(copied_network, resolution, weight)

    while community_size_combination(
            len(community) for community in communities) > community_size:
        subdivision = False
        for i, subdivided_community in enumerate(
                algorithm(copied_network.subgraph(communities[j]), resolution,
                          weight) for j in range(len(communities))):
            if len(subdivided_community) > 1:
                subdivision = True
                communities[i:i + 1] = subdivided_community

        if not subdivision:
            break

    return tuple(network.subgraph(community) for community in communities)


def get_measurements(
    network: nx.Graph,
    time: int,
    modification: str,
    site_combination: Optional[Callable[[Collection[float]], float]] = None
) -> tuple[float, ...]:
    """
    Returns the measurement distribution for a particular type of
    post-translational modification at a particular time of measurement.

    Args:
        network: The protein-protein interaction network.
        time: The time of measurement.
        modification: The type of post-translational modification.
        site_combination: An optional function to derive protein-specific
            measurements from site-specific measurements.

    Returns:
        Measurements for a particular type of post-translational
        modification at a particular time of measurement.
    """
    measurements = []
    for protein in network:
        sites = [
            network.nodes[protein][measurement]
            for measurement in network.nodes[protein]
            if len(measurement.split(" ")) == 3 and measurement.split(" ")[0] ==
            str(time) and measurement.split(" ")[1] == modification
        ]

        if sites:
            if site_combination:
                measurements.append(
                    math.log2(
                        site_combination(
                            [math.pow(2.0, site) for site in sites])))
            else:
                measurements.extend(sites)

    return tuple(measurements)


def get_measurement_enrichment(
    network: nx.Graph,
    communities: Iterable[nx.Graph],
    measurements: tuple[float, float] = (-1.0, 1.0),
    convert_measurement: Callable[
        [float, Collection[float]],
        float] = lambda measurement, measurements: measurement,
    site_combination: Optional[Callable[[Collection[float]], float]] = None,
    enrichment_test: Callable[
        [int, int, int, int],
        float] = lambda k, M, n, N: scipy.stats.hypergeom.sf(k - 1, M, n, N),
    multiple_testing_correction: Callable[
        [dict[tuple[nx.Graph, int, str], float]],
        dict[tuple[nx.Graph, int, str], float]] = correction.benjamini_hochberg,
) -> dict[nx.Graph, dict[int, dict[str, float]]]:
    """
    Test communities for enrichment of large protein-specific measurements for
    each time of measurement and type of post-translational modification with
    respect to the protein-protein interaction network.

    Args:
        network: The protein-protein interaction network.
        communities: The communities of the protein-protein interaction network.
        measurements: Proteins are classified by whether their representative
            measurement exceeds this range.
        convert_measurement: The function used convert the bounds to log2-fold
            measurements.
        site_combination: An optional function to derive protein-specific
            measurements from site-specific measurements.
        enrichment_test: The statistical test used to assess enrichment of
            proteins with large measurements by a community.
        multiple_testing_correction: The procedure to correct for multiple
            testing of multiple communities, times of measurement and types of
            post-translational modification.

    Returns:
        Corrected p-values for the enrichment of large protein-specific
        measurements by each community for each time of measurement and type of
        post-translational modification.
    """
    p_values = {}
    for time in get_times(network):
        for modification in get_modifications(network, time):
            measurement_range = (convert_measurement(
                measurements[0],
                get_measurements(network, time, modification,
                                 site_combination)),
                                 convert_measurement(
                                     measurements[1],
                                     get_measurements(network, time,
                                                      modification,
                                                      site_combination)))

            modified_proteins = len([
                measurement for measurement in get_measurements(
                    network, time, modification, site_combination)
                if measurement
            ])

            modified_community_proteins = [
                len([
                    measurement
                    for measurement in get_measurements(
                        community, time, modification, site_combination)
                    if measurement
                ])
                for community in communities
            ]

            target_proteins = len([
                measurement for measurement in get_measurements(
                    network, time, modification, site_combination)
                if measurement <= measurement_range[0] or
                measurement >= measurement_range[1]
            ])

            target_community_proteins = [
                len([
                    measurement
                    for measurement in get_measurements(
                        community, time, modification, site_combination)
                    if measurement <= measurement_range[0] or
                    measurement >= measurement_range[1]
                ])
                for community in communities
            ]

            p_values.update({(community, time, modification):
                             enrichment_test(target_community_proteins[i],
                                             modified_proteins, target_proteins,
                                             modified_community_proteins[i])
                             for i, community in enumerate(communities)})

    p_values = multiple_testing_correction(p_values)

    return {
        community: {
            time: {
                modification: p_values[(community, time, modification)]
                for modification in get_modifications(network, time)
            } for time in get_times(network)
        } for community in communities
    }


def get_measurement_location(
    network: nx.Graph,
    communities: Iterable[nx.Graph],
    site_combination: Optional[Callable[[Collection[float]], float]] = None,
    location_test: Callable[
        [Collection[float], Collection[float]],
        float] = lambda x, y: scipy.stats.ranksums(x, y).pvalue,
    multiple_testing_correction: Callable[
        [dict[tuple[nx.Graph, int, str], float]],
        dict[tuple[nx.Graph, int, str], float]] = correction.benjamini_hochberg,
) -> dict[nx.Graph, dict[int, dict[str, float]]]:
    """
    Test communities for difference tendencies in protein-specific measurements
    for each time of measurement and type of post-translational modification
    with respect to the remaining protein-protein interaction network.

    Args:
        network: The protein-protein interaction network.
        communities: The communities of the protein-protein interaction network.
        site_combination: An optional function to derive protein-specific
            measurements from site-specific measurements.
        location_test: The statistical test used to assess location of
            proteins with large measurements by a community.
        multiple_testing_correction: The procedure to correct for multiple
            testing of multiple communities, times of measurement and types of
            post-translational modification.

    Returns:
        Corrected p-values for the difference in central tendencies in
        protein-specific measurements of communities for each time of
        measurement and type of post-translational modification.
    """
    p_values = {}
    for time in get_times(network):
        for modification in get_modifications(network, time):
            network_measurements = [
                get_measurements(
                    nx.union_all([m
                                  for m in communities
                                  if m != community]), time, modification,
                    site_combination)
                for community in communities
            ]

            community_measurements = [
                get_measurements(community, time, modification,
                                 site_combination) for community in communities
            ]

            p_values.update({(community, time, modification):
                             location_test(community_measurements[i],
                                           network_measurements[i])
                             for i, community in enumerate(communities)
                             if community_measurements[i]})

    p_values = multiple_testing_correction(p_values)

    return {
        community: {
            time: {
                modification: p_values[(community, time, modification)]
                for modification in get_modifications(network, time)
                if (community, time, modification) in p_values
            } for time in get_times(network)
        } for community in communities
    }


def export(network: nx.Graph, basename: str) -> None:
    """
    Exports the protein-protein interaction network.

    Args:
        network: The protein-protein interaction network.
        basename: The base file name.
    """
    nx.write_graphml_xml(network,
                         f"{basename}.graphml",
                         named_key_ids=True,
                         infer_numeric_types=True)
