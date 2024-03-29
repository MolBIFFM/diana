"""protein-protein interaction network"""
import math
import os
import re
import statistics
from typing import (Callable, Collection, Container, Hashable, Iterable,
                    Literal, MutableSequence, Optional, Union)

import networkx as nx
import pandas as pd
import scipy.stats
from algorithms import correction, modularization
from databases import biogrid, corum, intact, mint, reactome, string, uniprot


def get_network() -> nx.Graph:
    """
    Returns an initialized protein-protein interaction network.

    Returns
        The protein-protein interaction network.

    """
    # Return an initialized protein-protein interaction network.
    return nx.Graph()


def add_proteins_from_table(
        network: nx.Graph,
        file: str,
        protein_accession_column: int | str,
        replicate_columns: Collection[int] | Collection[str],
        protein_accession_format: re.Pattern[str] = re.compile("^(.+)$"),
        replicate_format: re.Pattern[str] = re.compile("^(.+)$"),
        sheet_name: int | str = 0,
        header: int = 0,
        time: int = 0,
        modification: str = "PTM",
        number_replicates: int = 1,
        measurement_score: Callable[[float], float] = math.log2) -> None:
    """
    Parse UniProt protein accessions and protein-specific measurements from a
    tabular file and add proteins to a protein-protein interaction network.

    Args:
        network: The protein-protein interaction network.
        file: The file location of the file.
        protein_accession_column: The column containing UniProt protein
            accessions.
        replicate_columns: The columns containing replicates of measurements.
        protein_accession_format: A regular expression to extract protein
            accessions from a corresponding entry.
        replicate_format: A regular expression to extract measurements
            from a corresponding entry.
        sheet_name: The sheet to parse protein accessions from.
        header: The index of the header row.
        time: The time of measurement to associate with measurements.
        modification: An identifier for the type of post-translational
            modification to associate with measurements.
        number_replicates: The minimum number of replicates to accept a
            measurement.
        measurement_score: The function to convert the measurements
            reported to their binary logarithm.
    """
    # Parse proteins from a spreadsheet.
    if os.path.splitext(file)[1].lstrip(".") in (
            "xls",
            "xlsx",
            "xlsm",
            "xlsb",
            "odf",
            "ods",
            "odt",
    ):
        table = pd.read_excel(
            file,
            sheet_name=sheet_name,
            header=header,
            usecols=[protein_accession_column] + list(replicate_columns),
            dtype={
                protein_accession_column: str,
                **{column: str for column in replicate_columns},
            },
        )

    # Parse proteins from a tabular file.
    else:
        table = pd.read_table(
            file,
            header=header,
            usecols=[protein_accession_column] + list(replicate_columns),
            dtype={
                protein_accession_column: str,
                **{column: str for column in replicate_columns},
            },
        )

    # Parse protein accessions and measurements from the table asserting an
    # applicable replicate threshold.
    proteins: dict[str, tuple[float, ...]] = {}
    for _, row in table.iterrows():
        if pd.isna(row[protein_accession_column]):
            continue

        protein_accessions = [
            str(protein_accession) for protein_accession in
            protein_accession_format.findall(row[protein_accession_column])
        ]

        measurements = [
            float(replicate) for replicate_column in replicate_columns
            if not pd.isna(row[replicate_column])
            for replicate in replicate_format.findall(row[replicate_column])
        ]

        if len(measurements) >= min(number_replicates, len(replicate_columns)):
            for protein_accession in protein_accessions:
                # https://www.uniprot.org/help/accession_numbers
                # https://www.uniprot.org/help/alternative_products
                if re.fullmatch(
                        r"([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]"
                        r"([A-Z][A-Z0-9]{2}[0-9]){1,2})(-[1-9][0-9]+)?",
                        protein_accession):
                    proteins[protein_accession] = tuple(
                        measurement_score(measurement)
                        for measurement in measurements)

    # Add protein accessions to the protein-protein interaction network and
    # annotate nodes with measurements.
    network.add_nodes_from(proteins)
    for protein, replicates in proteins.items():
        for r, replicate in enumerate(replicates, start=1):
            network.nodes[protein][f"{time} {modification} R{r}"] = replicate


def add_sites_from_table(
    network: nx.Graph,
    file: str,
    protein_accession_column: int | str,
    position_column: int | str,
    replicate_columns: Collection[int] | Collection[str],
    protein_accession_format: re.Pattern[str] = re.compile("^(.+)$"),
    position_format: re.Pattern[str] = re.compile("^(.+)$"),
    replicate_format: re.Pattern[str] = re.compile("^(.+)$"),
    sheet_name: int | str = 0,
    header: int = 0,
    time: int = 0,
    modification: str = "PTM",
    number_sites: int = 100,
    number_replicates: int = 1,
    replicate_average: Callable[[Iterable[float]], float] = statistics.mean,
    measurement_score: Callable[[float], float] = math.log2,
    site_prioritization: Callable[[float],
                                  float] = lambda site: abs(math.log2(site)),
    site_order: Callable[[tuple[int, float]],
                         Union[int, float]] = lambda site: site[0]
) -> None:
    """
    Parse UniProt protein accessions and site-specific measurements from a
    tabular file and add proteins to a protein-protein interaction network.

    Args:
        network: The protein-protein interaction network.
        file: The file location of the file.
        protein_accession_column: The column containing UniProt protein
            accessions.
        position_column: The column containing sites corresponding to
            measurements.
        replicate_columns: The columns containing replicates of measurements.
        protein_accession_format: A regular expression to extract protein
            accessions from a corresponding entry.
        position_format: A regular expression to extract measurement positions
            such as modification sites from a corresponding entry.
        replicate_format: A regular expression to extract measurements
            from a corresponding entry.
        sheet_name: The sheet to parse protein accessions from.
        header: The index of the header row.
        time: The time of measurement to associate with measurements.
        modification: An identifier for the type of post-translational
            modification to associate with measurements.
        number_sites: The maximum number of measurements to associate with a
            protein-accession, prioritized by value.
        number_replicates: The minimum number of replicates to accept a
            measurement.
        replicate_average: A function to derive site-specific measurements from
            replicates for site prioritization.
        measurement_score: The function to convert the measurements
            reported to their binary logarithm.
        site_prioritization: A function of a measurement to derive a
            representative value for site prioritization.
        site_order: A function of a measurement to derive a
            representative value for site order.
    """
    # Parse proteins from a spreadsheet.
    if os.path.splitext(file)[1].lstrip(".") in (
            "xls",
            "xlsx",
            "xlsm",
            "xlsb",
            "odf",
            "ods",
            "odt",
    ):
        table = pd.read_excel(
            file,
            sheet_name=sheet_name,
            header=header,
            usecols=[protein_accession_column, position_column] +
            list(replicate_columns),
            dtype={
                protein_accession_column: str,
                position_column: str,
                **{column: str for column in replicate_columns},
            },
        )

    # Parse proteins from a tabular file.
    else:
        table = pd.read_table(
            file,
            header=header,
            usecols=[protein_accession_column, position_column] +
            list(replicate_columns),
            dtype={
                protein_accession_column: str,
                position_column: str,
                **{column: str for column in replicate_columns},
            },
        )

    # Parse protein accessions, positions of modification sites, and
    # measurements from the table asserting an applicable replicate threshold.
    proteins: dict[str, dict[int, tuple[float, ...]]] = {}
    for _, row in table.iterrows():
        if pd.isna(row[protein_accession_column]):
            continue

        protein_accessions = [
            str(protein_accession) for protein_accession in
            protein_accession_format.findall(row[protein_accession_column])
        ]

        if not pd.isna(row[position_column]):
            positions = [
                int(position)
                for position in position_format.findall(row[position_column])
            ]
        else:
            positions = []

        # The row contains more protein accessions than positions of
        # modification sites.
        if len(protein_accessions) > len(positions):
            positions.extend(
                [0 for _ in range(len(positions), len(protein_accessions))])

        # The row contains less protein accessions than positions of
        # modification sites.
        elif len(protein_accessions) < len(positions):
            positions = positions[:len(protein_accessions)]

        measurements = [
            float(replicate) for replicate_column in replicate_columns
            if not pd.isna(row[replicate_column])
            for replicate in replicate_format.findall(row[replicate_column])
        ]

        if len(measurements) >= min(number_replicates, len(replicate_columns)):
            for protein_accession, position in zip(protein_accessions,
                                                   positions):
                # https://www.uniprot.org/help/accession_numbers
                # https://www.uniprot.org/help/alternative_products
                if re.fullmatch(
                        r"([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]"
                        r"([A-Z][A-Z0-9]{2}[0-9]){1,2})(-[1-9][0-9]+)?",
                        protein_accession):
                    if protein_accession not in proteins:
                        proteins[protein_accession] = {}

                    proteins[protein_accession][position] = tuple(
                        measurement_score(measurement)
                        for measurement in measurements)

    # Add protein accessions to the protein-protein interaction network and
    # annotate nodes with measurements.
    network.add_nodes_from(proteins)
    for protein, sites in proteins.items():
        for s, (_, replicates) in enumerate(sorted(
                sorted(
                    sites.items(),
                    key=lambda site: site_prioritization(
                        replicate_average(
                            [math.pow(2.0, replicate)
                             for replicate in site[1]])))[-number_sites:],
                key=lambda site: site_order(
                    (site[0],
                     math.pow(
                         2.0,
                         replicate_average([
                             math.pow(2.0, replicate) for replicate in site[1]
                         ]))))),
                                            start=1):
            for m, measurement in enumerate(replicates, start=1):
                network.nodes[protein][
                    f"{time} {modification} S{s} R{m}"] = measurement


def map_proteins(network: nx.Graph,
                 organism: int = 9606,
                 file: Optional[str] = None) -> frozenset[str]:
    """
    Map proteins in a protein-protein interaction network to their primary
    UniProt identifiers and remove proteins not present in Swiss-Prot for the
    organism of interest. Isoform identifiers are maintained on primary, but not
    transferred to secondary UniProt accessions.

    Args:
        network: The protein-protein interaction network to map proteins from.
        organism: The NCBI taxonomy identifier for the organism of interest.
        file: The optional local file location to parse accessions from.

    Returns:
        The proteins not present in the queried portion of Swiss-Prot specific
        to the organism of interest.
    """
    # Compile a map from UniProt protein accessions present in the
    # protein-protein interaction network to primary UniProt protein accessions
    # and associated gene and protein names.
    mapping, gene_name, protein_name = {}, {}, {}
    for accessions, gene, protein in uniprot.get_swiss_prot_entries(
            organism, file):
        for node in network:
            if node.split("-")[0] in accessions:
                if "-" in node and accessions.index(node.split("-")[0]) == 0:
                    mapping[node] = node
                    gene_name[node] = gene
                    protein_name[node] = protein
                else:
                    mapping[node] = accessions[0]
                    gene_name[accessions[0]] = gene
                    protein_name[accessions[0]] = protein

    # Relabel nodes of the protein-protein interaction network to their
    # associated primary UniProt protein accession.
    nx.relabel_nodes(network, mapping, copy=False)

    # Annotate nodes of the protein-protein interaction network with
    for protein, name in protein_name.items():
        if protein in network:
            network.nodes[protein]["protein"] = name
            network.nodes[protein]["gene"] = gene_name[protein]

    # Return nodes not associated with primary UniProt protein accessions.
    return frozenset(set(network).difference(mapping))


def get_proteins(
    network: nx.Graph,
    time: int,
    modification: str,
    site_average: Optional[Callable[[Iterable[float]], float]] = lambda sites:
    max(sites, key=lambda site: abs(math.log2(site))),
    replicate_average: Optional[Callable[[Iterable[float]],
                                         float]] = statistics.mean,
    measurement_range: tuple[float, float] = (0.0, 0.0)
) -> frozenset[str]:
    """
    Returns proteins of a protein-protein interaction network exceeding a
    specified range of measurements.

    Args:
        network: The protein-protein interaction network.
        time: The time of measurement.
        modification: The type of post-translational modification.
        site_average: A function to derive protein-specific measurements from
            site-specific measurements.
        replicate_average: A function to derive site-specific measurements from
            replicates.
        measurement_range: A range that is to be exceeded by average
            measurements of the proteins.

    Returns:
        Proteins the average measurement of which for a particular type of
        post-translational modification at a particular time of measurement
        exceeds the specified range.
    """
    # Return the set of proteins from the protein-protein interaction network
    # the average site or replicate measurements of which for a particular type
    # of post-translational modification and time of measurement are all at most
    # a lower or at least an upper threshold.
    proteins = []
    for protein in network:
        sites = [[
            network.nodes[protein][attribute]
            for attribute in network.nodes[protein]
            if re.fullmatch(fr"{time} {modification} (S{s+1} )?R[0-9]+",
                            attribute)
        ]
                 for s in range(get_sites(network, time, modification, protein))
                ]

        if site_average is not None and replicate_average is not None:
            if sites and not (measurement_range[0] < math.log2(
                    site_average([
                        replicate_average(
                            [math.pow(2.0, replicate)
                             for replicate in site])
                        for site in sites
                    ])) < measurement_range[1]):
                proteins.append(protein)

        elif replicate_average is not None:
            if sites and not all(measurement_range[0] < math.log2(
                    replicate_average(
                        [math.pow(2.0, replicate)
                         for replicate in site])) < measurement_range[1]
                                 for site in sites):
                proteins.append(protein)
        else:
            if sites and not all(
                    measurement_range[0] < replicate < measurement_range[1]
                    for site in sites
                    for replicate in site):
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
    # Return times of measurement extracted from identifiers for measurements
    # from node annotations of a protein-protein interaction network.
    return tuple(
        sorted(
            set(
                int(attribute.split(" ")[0])
                for protein in network
                for attribute in network.nodes[protein]
                if re.fullmatch(r"[0-9]+ [A-Za-z]+ (S[0-9]+ )?R[0-9]+",
                                attribute))))


def get_modifications(network: nx.Graph,
                      time: int,
                      proteins: bool = True,
                      sites: bool = True) -> tuple[str, ...]:
    """
    Returns the types of post-translational modification represented in a
    protein-protein interaction network at a particular time of measurement.

    Args:
        network: The protein-protein interaction network.
        time: The time of measurement.
        proteins: If true, consider protein-specific types of post-translational
            modification.
        sites: If true, consider site-specific types of post-translational
            modification.

    Returns:
        The types of post-translational modification associated with any
        measurements of any protein in the protein-protein interaction network
        at a particular time of measurement.
    """
    # Each modification is either protein- or site-specific.
    if not (proteins or sites):
        return ()

    # Return identifiers for types of post-translational modification extracted
    # from identifiers for measurements from node annotations of a
    # protein-protein interaction network.
    return tuple(
        sorted(
            set(
                attribute.split(" ")[1]
                for protein in network
                for attribute in network.nodes[protein]
                if re.fullmatch(
                    fr"{time} [A-Za-z]+ (S[0-9]+ )?R[0-9]+"
                    if proteins and sites else (
                        fr"{time} [A-Za-z]+ R[0-9]+" if proteins else
                        fr"{time} [A-Za-z]+ S[0-9]+ R[0-9]+"), attribute))))


def is_modification(network: nx.Graph,
                    time: int,
                    modification: Optional[str],
                    proteins: bool = True,
                    sites: bool = True) -> bool:
    """
    Indicate whether a type of post-translational modification is represented in
    a protein-protein interaction network at a particular time of measurement.

    Args:
        network: The protein-protein interaction network.
        time: The time of measurement.
        modification: The type of post-translational modification.
        proteins: If true, consider protein-specific types of post-translational
            modification.
        sites: If true, consider site-specific types of post-translational
            modification.

    Returns:
        True if the type of post-translational modification is represented in
        the protein-protein interaction network at the time of measurement, else
        false.
    """
    # Each modification is either protein- or site-specific.
    if not (proteins or sites):
        return False

    # Indicate whether the identifier corresponds to an identifier of a protein-
    # or site-specific type of post-translational modification represented in
    # the protein-protein interaction network.
    return any(
        re.fullmatch(
            fr"{time} {modification} (S[0-9]+ )?R[0-9]+" if proteins and sites
            else (fr"{time} {modification} R[0-9]+" if proteins else
                  fr"{time} [A-Za-z]+ S[0-9]+ R[0-9]+"), attribute)
        for protein in network
        for attribute in network.nodes[protein])


def get_sites(network: nx.Graph, time: int, modification: str,
              protein: str) -> int:
    """
    Returns the number of sites of a protein represented in a protein-protein
    interaction network for a particular type of post-translational modification
    at a particular time of measurement for site-specific modifications and 1
    for protein-specific modifications.

    Args:
        network: The protein-protein interaction network.
        time: The time of measurement.
        modification: The type of post-translational modification.
        protein: The protein.

    Returns:
        The number of sites of the protein in a protein-protein interaction
        network for a particular type of post-translational modification at a
        particular time of measurement.
    """
    # Return the maximum number of modification sites for a particular type of
    # post-translational modification and time of modification associated with
    # any protein represented in the protein-protein interaction network.
    return max(
        int(attribute.split(" ")[2].replace("S", "")) if re.
        fullmatch(fr"{time} {modification} S[0-9]+ R[0-9]+", attribute) else
        (1 if re.fullmatch(fr"{time} {modification} R[0-9]+", attribute) else 0)
        for attribute in network.nodes[protein])


def set_measurements(
    network: nx.Graph, site_average: dict[str, Callable[[Iterable[float]],
                                                        float]],
    replicate_average: dict[str, Callable[[Iterable[float]], float]],
    measurements: dict[str, tuple[float, float]],
    measurement_score: dict[str, Callable[[float, Collection[float]], float]]
) -> None:
    """
    Annotate proteins in a protein-protein interaction network with aggregates
    and categorization of measurements.

    Args:
        network: The protein-protein interaction network.
        site_average: Modification-specific functions to derive protein-specific
            measurements from site-specific measurements.
        replicate_average: Modification-specific functions to derive
            site-specific measurements from replicate measurements.
        measurements: Modification-specific range of averaged measurements by
            which proteins are categorized.
        measurement_score: Modification-specific functions used convert the
            measurement range bounds to their binary logarithm.
    """
    # Compile the set of times of measurement represented in the protein-protein
    # interaction network.
    times = get_times(network)

    # Compile a map from times of measurement to types of post-translational
    # modification represented in the protein-protein interaction network.
    modifications = {time: get_modifications(network, time) for time in times}

    # Compile a map from times of measurement and types of post-translational
    # modifications to measurement ranges.
    measurement_range = {
        time: {
            modification:
            (measurement_score.get(
                modification, lambda measurement, measurements: measurement)
             (
                 measurements.get(modification, (-1.0, 1.0))[0],
                 get_measurements(
                     network, time, modification,
                     site_average.get(
                         modification, lambda sites: max(
                             sites, key=lambda site: abs(math.log2(site)))),
                     replicate_average.get(modification, statistics.mean)),
             ),
             measurement_score.get(
                 modification,
                 lambda measurement, measurements: measurement)(
                     measurements.get(
                         modification,
                         (-1.0, 1.0))[1],
                     get_measurements(
                         network, time, modification,
                         site_average.get(
                             modification,
                             lambda sites: max(
                                 sites, key=lambda site: abs(math.log2(site)))),
                         replicate_average.get(modification, statistics.mean))))
            for modification in modifications[time]
        } for time in times
    }

    # Annotate nodes of the protein-protein interaction network with
    # categorization of
    for time in times:
        for prt in network:
            summary = {}
            for modification in modifications[time]:
                sites = [[
                    network.nodes[prt][attribute] for attribute in
                    network.nodes[prt] if re.fullmatch(
                        fr"{time} {modification} (S{s+1} )?R[0-9]+", attribute)
                ] for s in range(get_sites(network, time, modification, prt))]

                if sites:
                    if is_modification(network,
                                       time,
                                       modification,
                                       proteins=False):
                        for s, site in enumerate(sites, start=1):
                            network.nodes[prt][
                                f"{time} {modification} S{s}"] = math.log2(
                                    replicate_average.get(
                                        modification, statistics.mean)([
                                            math.pow(2.0, replicate)
                                            for replicate in site
                                        ]))

                    network.nodes[prt][f"{time} {modification}"] = math.log2(
                        site_average.get(
                            modification, lambda sites: max(
                                sites, key=lambda site: abs(math.log2(site))))([
                                    replicate_average.get(
                                        modification, statistics.mean)([
                                            math.pow(2.0, replicate)
                                            for replicate in site
                                        ])
                                    for site in sites
                                ]))

                    if (network.nodes[prt][f"{time} {modification}"] >= max(
                            0.0, measurement_range[time][modification][1])):
                        summary[modification] = "UP"
                    elif (max(0.0, measurement_range[time][modification][1]) >
                          network.nodes[prt][f"{time} {modification}"] >= max(
                              0.0, 0.5 *
                              (0.5 *
                               (measurement_range[time][modification][0] +
                                measurement_range[time][modification][1]) +
                               measurement_range[time][modification][1]))):
                        summary[modification] = "MID-UP"
                    elif (min(0.0, measurement_range[time][modification][0]) <
                          network.nodes[prt][f"{time} {modification}"] <= min(
                              0.0, 0.5 *
                              (measurement_range[time][modification][0] + 0.5 *
                               (measurement_range[time][modification][0] +
                                measurement_range[time][modification][1])))):
                        summary[modification] = "MID-DOWN"
                    elif (network.nodes[prt][f"{time} {modification}"] <= min(
                            0.0, measurement_range[time][modification][0])):
                        summary[modification] = "DOWN"
                    else:
                        summary[modification] = "MID"

            if summary:
                network.nodes[prt][str(time)] = " ".join(
                    f"{modification} {category}"
                    for modification, category in summary.items())
            else:
                network.nodes[prt][str(time)] = ""


def get_neighbors_from_biogrid(
        network: nx.Graph,
        experimental_system: Optional[Container[str]] = None,
        experimental_system_type: Optional[Container[str]] = None,
        interaction_throughput: Optional[Container[str]] = None,
        multi_validated_physical: bool = False,
        organism: int = 9606,
        version: Optional[str] = None,
        file: Optional[str] = None,
        file_uniprot: Optional[str] = None) -> set[str]:
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
        file: The optional local file location to parse interactions from.
        file_uniprot: The optional local file location to parse accessions from.

    Returns:
        Neighbors of the protein-protein interaction network in BioGRID.
    """
    # Return the set of proteins from BioGRID not represented but interacting
    # with proteins in the protein-protein interaction network.
    neighbors = set()
    for interactor_a, interactor_b in biogrid.get_protein_interactions(
            experimental_system, experimental_system_type,
            interaction_throughput, multi_validated_physical, organism, version,
            file, file_uniprot):
        # https://www.uniprot.org/help/accession_numbers
        # https://www.uniprot.org/help/alternative_products
        if (interactor_a in network and interactor_b not in network and
                re.fullmatch(
                    r"([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]"
                    r"([A-Z][A-Z0-9]{2}[0-9]){1,2})(-[1-9][0-9]+)?",
                    interactor_b)):
            neighbors.add(interactor_b)

        elif (interactor_a not in network and interactor_b in network and
              re.fullmatch(
                  r"([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]"
                  r"([A-Z][A-Z0-9]{2}[0-9]){1,2})(-[1-9][0-9]+)?",
                  interactor_a)):
            neighbors.add(interactor_a)
    return neighbors


def add_protein_interactions_from_biogrid(
        network: nx.Graph,
        experimental_system: Optional[Container[str]] = None,
        experimental_system_type: Optional[Container[str]] = None,
        interaction_throughput: Optional[Container[str]] = None,
        multi_validated_physical: bool = False,
        organism: int = 9606,
        version: Optional[str] = None,
        file: Optional[str] = None,
        file_uniprot: Optional[str] = None) -> None:
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
        file: The optional local file location to parse interactions from.
        file_uniprot: The optional local file location to parse accessions from.
    """
    # Add protein-protein interactions between different proteins from BioGRID
    # to the protein-protein interaction network.
    for interactor_a, interactor_b in biogrid.get_protein_interactions(
            experimental_system, experimental_system_type,
            interaction_throughput, multi_validated_physical, organism, version,
            file, file_uniprot):
        if (interactor_a in network and interactor_b in network and
                interactor_a != interactor_b):
            network.add_edge(interactor_a, interactor_b)
            network.edges[interactor_a, interactor_b]["BioGRID"] = 1.0


def get_neighbors_from_corum(network: nx.Graph,
                             purification_methods: Optional[
                                 Container[str]] = None,
                             organism: int = 9606,
                             file: Optional[str] = None,
                             file_uniprot: Optional[str] = None) -> set[str]:
    """
    Returns proteins interacting with proteins in a protein-protein interaction
    network from CORUM.

    Args:
        network: The protein-protein interaction network.
        purification_methods: The accepted PSI-MI identifiers or terms for the
            protein complex purification method. If none are specified, any are
            accepted.
        organism: The NCBI taxonomy identifier for the organism of interest.
        file: The optional local file location to parse interactions from.
        file_uniprot: The optional local file location to parse accessions from.

    Returns:
        Neighbors of the protein-protein interaction network in CORUM.
    """
    # Return the set of proteins from CORUM not represented but interacting with
    # proteins in the protein-protein interaction network.
    neighbors = set()
    for interactor_a, interactor_b in corum.get_protein_interactions(
            purification_methods, organism, file, file_uniprot):
        # https://www.uniprot.org/help/accession_numbers
        # https://www.uniprot.org/help/alternative_products
        if (interactor_a in network and interactor_b not in network and
                re.fullmatch(
                    r"([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]"
                    r"([A-Z][A-Z0-9]{2}[0-9]){1,2})(-[1-9][0-9]+)?",
                    interactor_b)):
            neighbors.add(interactor_b)

        elif (interactor_a not in network and interactor_b in network and
              re.fullmatch(
                  r"([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]"
                  r"([A-Z][A-Z0-9]{2}[0-9]){1,2})(-[1-9][0-9]+)?",
                  interactor_a)):
            neighbors.add(interactor_a)

    return neighbors


def add_protein_interactions_from_corum(
        network: nx.Graph,
        purification_methods: Optional[Container[str]] = None,
        organism: int = 9606,
        file: Optional[str] = None,
        file_uniprot: Optional[str] = None) -> None:
    """
    Adds protein-protein interactions from CORUM to a protein-protein
    interaction network.

    Args:
        network: The protein-protein interaction network.
        purification_methods: The accepted PSI-MI identifiers or terms for the
            protein complex purification method. If none are specified, any are
            accepted.
        organism: The NCBI taxonomy identifier for the organism of interest.
        file: The optional local file location to parse interactions from.
        file_uniprot: The optional local file location to parse accessions from.
    """
    # Add protein-protein interactions between different proteins from CORUM to
    # the protein-protein interaction network.
    for interactor_a, interactor_b in corum.get_protein_interactions(
            purification_methods, organism, file, file_uniprot):
        if (interactor_a in network and interactor_b in network and
                interactor_a != interactor_b):
            network.add_edge(interactor_a, interactor_b)
            network.edges[interactor_a, interactor_b]["CORUM"] = 1.0


def get_neighbors_from_intact(
        network: nx.Graph,
        interaction_detection_methods: Optional[Container[str]] = None,
        interaction_types: Optional[Container[str]] = None,
        psi_mi_score: float = 0.0,
        organism: int = 9606,
        file: Optional[str] = None,
        file_uniprot: Optional[str] = None) -> set[str]:
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
        file: The optional local file location to parse interactions from.
        file_uniprot: The optional local file location to parse accessions from.

    Returns:
        Neighbors of the protein-protein interaction network in IntAct.
    """
    # Return the set of proteins from IntAct not represented but interacting
    # with proteins in the protein-protein interaction network.
    neighbors = set()
    for interactor_a, interactor_b, _ in intact.get_protein_interactions(
            interaction_detection_methods, interaction_types, psi_mi_score,
            organism, file, file_uniprot):
        # https://www.uniprot.org/help/accession_numbers
        # https://www.uniprot.org/help/alternative_products
        if (interactor_a in network and interactor_b not in network and
                re.fullmatch(
                    r"([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]"
                    r"([A-Z][A-Z0-9]{2}[0-9]){1,2})(-[1-9][0-9]+)?",
                    interactor_b)):
            neighbors.add(interactor_b)

        elif (interactor_a not in network and interactor_b in network and
              re.fullmatch(
                  r"([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]"
                  r"([A-Z][A-Z0-9]{2}[0-9]){1,2})(-[1-9][0-9]+)?",
                  interactor_a)):
            neighbors.add(interactor_a)

    return neighbors


def add_protein_interactions_from_intact(
        network: nx.Graph,
        interaction_detection_methods: Optional[Container[str]] = None,
        interaction_types: Optional[Container[str]] = None,
        psi_mi_score: float = 0.0,
        organism: int = 9606,
        file: Optional[str] = None,
        file_uniprot: Optional[str] = None) -> None:
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
        file: The optional local file location to parse interactions from.
        file_uniprot: The optional local file location to parse accessions from.
    """
    # Add protein-protein interactions between different proteins from IntAct to
    # the protein-protein interaction network.
    for interactor_a, interactor_b, score in intact.get_protein_interactions(
            interaction_detection_methods, interaction_types, psi_mi_score,
            organism, file, file_uniprot):
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
                                Container[str]] = None,
                            interaction_types: Optional[Container[str]] = None,
                            psi_mi_score: float = 0.0,
                            organism: int = 9606,
                            file: Optional[str] = None,
                            file_uniprot: Optional[str] = None) -> set[str]:
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
        file: The optional local file location to parse interactions from.
        file_uniprot: The optional local file location to parse accessions from.

    Returns:
        Neighbors of the protein-protein interaction network in MINT.
    """
    # Return the set of proteins from MINT not represented but interacting with
    # proteins in the protein-protein interaction network.
    neighbors = set()
    for interactor_a, interactor_b, _ in mint.get_protein_interactions(
            interaction_detection_methods, interaction_types, psi_mi_score,
            organism, file, file_uniprot):
        # https://www.uniprot.org/help/accession_numbers
        # https://www.uniprot.org/help/alternative_products
        if (interactor_a in network and interactor_b not in network and
                re.fullmatch(
                    r"([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]"
                    r"([A-Z][A-Z0-9]{2}[0-9]){1,2})(-[1-9][0-9]+)?",
                    interactor_b)):
            neighbors.add(interactor_b)

        elif (interactor_a not in network and interactor_b in network and
              re.fullmatch(
                  r"([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]"
                  r"([A-Z][A-Z0-9]{2}[0-9]){1,2})(-[1-9][0-9]+)?",
                  interactor_a)):
            neighbors.add(interactor_a)

    return neighbors


def add_protein_interactions_from_mint(
        network: nx.Graph,
        interaction_detection_methods: Optional[Container[str]] = None,
        interaction_types: Optional[Container[str]] = None,
        psi_mi_score: float = 0.0,
        organism: int = 9606,
        file: Optional[str] = None,
        file_uniprot: Optional[str] = None) -> None:
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
        file: The optional local file location to parse interactions from.
        file_uniprot: The optional local file location to parse accessions from.
    """
    # Add protein-protein interactions between different proteins from MINT to
    # the protein-protein interaction network.
    for interactor_a, interactor_b, score in mint.get_protein_interactions(
            interaction_detection_methods, interaction_types, psi_mi_score,
            organism, file, file_uniprot):
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
        organism: int = 9606,
        file: Optional[str] = None,
        file_uniprot: Optional[str] = None) -> set[str]:
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
        file: The optional local file location to parse interactions from.
        file_uniprot: The optional local file location to parse accessions from.

    Returns:
        Neighbors of the protein-protein interaction network in Reactome.
    """
    # Return the set of proteins from Reactome not represented but interacting
    # with proteins in the protein-protein interaction network.
    neighbors = set()
    for interactor_a, interactor_b in reactome.get_protein_interactions(
            interaction_type, interaction_context, organism, file,
            file_uniprot):
        # https://www.uniprot.org/help/accession_numbers
        # https://www.uniprot.org/help/alternative_products
        if (interactor_a in network and interactor_b not in network and
                re.fullmatch(
                    r"([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]"
                    r"([A-Z][A-Z0-9]{2}[0-9]){1,2})(-[1-9][0-9]+)?",
                    interactor_b)):
            neighbors.add(interactor_b)

        elif (interactor_a not in network and interactor_b in network and
              re.fullmatch(
                  r"([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]"
                  r"([A-Z][A-Z0-9]{2}[0-9]){1,2})(-[1-9][0-9]+)?",
                  interactor_a)):
            neighbors.add(interactor_a)

    return neighbors


def add_protein_interactions_from_reactome(
        network: nx.Graph,
        interaction_type: Optional[Container[str]] = None,
        interaction_context: Optional[Container[str]] = None,
        organism: int = 9606,
        file: Optional[str] = None,
        file_uniprot: Optional[str] = None) -> None:
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
        file: The optional local file location to parse interactions from.
        file_uniprot: The optional local file location to parse accessions from.
    """
    # Add protein-protein interactions between different proteins from Reactome
    # to the protein-protein interaction network.
    for interactor_a, interactor_b in reactome.get_protein_interactions(
            interaction_type, interaction_context, organism, file,
            file_uniprot):
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
                              version: float = 11.5,
                              any_score: bool = False,
                              file: Optional[str] = None,
                              file_accession_map: Optional[str] = None,
                              file_uniprot: Optional[str] = None) -> set[str]:
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
        any_score: If true, include a protein-protein interaction if it meets
            any threshold else only if it meets all.
        file: The optional local file location to parse interactions from.
        file_accession_map: The optional local file location to parse UniProt
            accessions from.
        file_uniprot: The optional local file location to parse accessions from.

    Returns:
        Neighbors of the protein-protein interaction network in STRING.
    """
    # Return the set of proteins from STRING not represented but interacting
    # with proteins in the protein-protein interaction network.
    neighbors = set()
    for interactor_a, interactor_b, _ in string.get_protein_interactions(
            neighborhood, neighborhood_transferred, fusion, cooccurence,
            homology, coexpression, coexpression_transferred, experiments,
            experiments_transferred, database, database_transferred, textmining,
            textmining_transferred, combined_score, physical, organism, version,
            any_score, file, file_accession_map, file_uniprot):
        # https://www.uniprot.org/help/accession_numbers
        # https://www.uniprot.org/help/alternative_products
        if (interactor_a in network and interactor_b not in network and
                re.fullmatch(
                    r"([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]"
                    r"([A-Z][A-Z0-9]{2}[0-9]){1,2})(-[1-9][0-9]+)?",
                    interactor_b)):
            neighbors.add(interactor_b)

        elif (interactor_a not in network and interactor_b in network and
              re.fullmatch(
                  r"([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]"
                  r"([A-Z][A-Z0-9]{2}[0-9]){1,2})(-[1-9][0-9]+)?",
                  interactor_a)):
            neighbors.add(interactor_a)

    return neighbors


def add_protein_interactions_from_string(
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
        organism: int = 9606,
        version: float = 11.5,
        any_score: bool = False,
        file: Optional[str] = None,
        file_accession_map: Optional[str] = None,
        file_uniprot: Optional[str] = None) -> None:
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
        any_score: If true, include a protein-protein interaction if it meets
            any threshold else only if it meets all.
        file: The optional local file location to parse interactions from.
        file_accession_map: The optional local file location to parse UniProt
            accessions from.
        file_uniprot: The optional local file location to parse accessions from.
    """
    # Add protein-protein interactions between different proteins from STRING
    # to the protein-protein interaction network.
    for interactor_a, interactor_b, score in string.get_protein_interactions(
            neighborhood, neighborhood_transferred, fusion, cooccurence,
            homology, coexpression, coexpression_transferred, experiments,
            experiments_transferred, database, database_transferred, textmining,
            textmining_transferred, combined_score, physical, organism, version,
            any_score, file, file_accession_map, file_uniprot):
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


def get_databases(
    network: nx.Graph
) -> tuple[Literal["BioGRID", "CORUM", "IntAct", "MINT", "Reactome", "STRING"],
           ...]:
    """
    Returns the protein-protein interaction databases represented in a
    protein-protein interaction network.

    Args:
        network: The protein-protein interaction network.

    Returns:
        The databases of protein-protein interaction confidence scores
        associated with edges in the protein-protein interaction network.
    """
    # Return protein-protein interaction databases associated with edges of the
    # protein-protein interaction network.
    return tuple(
        sorted(
            set(database for edge in network.edges()
                for database in network.edges[edge]).intersection((
                    "BioGRID",
                    "CORUM",
                    "IntAct",
                    "MINT",
                    "Reactome",
                    "STRING",
                ))))


def set_edge_weights(
    network: nx.Graph,
    weight: Callable[[
        dict[Literal["BioGRID", "CORUM", "IntAct", "MINT", "Reactome",
                     "STRING"], float]
    ], float] = lambda confidence_scores: int(bool(confidence_scores.values())),
    attribute: str = "weight",
) -> None:
    """
    Set average edge weights of a protein-protein interaction network.

    Args:
        network: The protein-protein interaction network.
        weight: The function to derive a average confidence score from
            database-specific scores.
        attribute: The attribute name of the weight.
    """
    # Compile a set of protein-protein interaction databases associated with
    # edges of the protein-protein interaction network.
    databases = get_databases(network)

    # Annotate edges of the protein-protein interaction network with average
    # confidence
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
    # Remove attribute from all edges of the protein-protein interaction
    # network.
    for _, _, data in network.edges(data=True):
        if attribute in data:
            del data[attribute]


def get_communities(
        network: nx.Graph,
        community_size: int,
        community_size_average: Callable[[Iterable[int]],
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
        community_size_average: The function to derive a decisive value to from
            meeting the threshold from the community sizes.
        algorithm: The community detection algorithm.
        resolution: The resolution parameter for modularity.
        weight: The attribute name of the edge weight.

    Returns:
        Modules of the protein-protein interaction network.
    """
    # No communities can be detected as the network contains no edges.
    if not network.number_of_edges():
        return tuple()

    # Identify communities disregarding isolate nodes.
    copied_network = network.copy()
    copied_network.remove_nodes_from(list(nx.isolates(copied_network)))

    # Apply the community detection algorithm.
    communities = algorithm(copied_network, resolution, weight)

    # Iteratively subdivide communities until average module size threshold is
    # met.
    while community_size_average(
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

    # Return communities as induced subgraphs of the protein-protein interaction
    # network.
    return tuple(network.subgraph(community) for community in communities)


def get_measurements(
    network: nx.Graph,
    time: int,
    modification: str,
    site_average: Optional[Callable[[Iterable[float]], float]] = None,
    replicate_average: Optional[Callable[[Iterable[float]], float]] = None
) -> tuple[float, ...]:
    """
    Returns the measurement distribution for a particular type of
    post-translational modification at a particular time of measurement.

    Args:
        network: The protein-protein interaction network.
        time: The time of measurement.
        modification: The type of post-translational modification.
        site_average: An optional function to derive protein-specific
            measurements from site-specific measurements.
        replicate_average: The function to derive site-specific measurements
            from replicates.

    Returns:
        Measurements for a particular type of post-translational
        modification at a particular time of measurement.
    """
    # Compile the distribution of the average measurement of a particular type
    # of post-translational modification at a particular time of measurement.
    measurements = []
    for protein in network:
        sites = [[
            network.nodes[protein][attribute]
            for attribute in network.nodes[protein]
            if re.fullmatch(fr"{time} {modification} (S{s+1} )?R[0-9]+",
                            attribute)
        ]
                 for s in range(get_sites(network, time, modification, protein))
                ]

        if sites:
            if site_average is not None and replicate_average is not None:
                measurements.append(
                    math.log2(
                        site_average([
                            replicate_average([
                                math.pow(2.0, replicate) for replicate in site
                            ]) for site in sites
                        ])))
            elif replicate_average is not None:
                measurements.extend([
                    math.log2(
                        replicate_average(
                            [math.pow(2.0, replicate)
                             for replicate in site]))
                    for site in sites
                ])
            else:
                for site in sites:
                    measurements.extend(site)

    return tuple(measurements)


def get_enrichment(
    network: nx.Graph,
    communities: Iterable[nx.Graph],
    measurement_ranges: dict[str, tuple[float, float]],
    measurement_score: dict[str, Callable[[float, Collection[float]], float]],
    site_average: dict[str, Optional[Callable[[Iterable[float]], float]]],
    replicate_average: dict[str, Optional[Callable[[Iterable[float]], float]]],
    enrichment_test: Callable[[int, int, int, int],
                              float] = lambda k, M, n, N: float(
                                  scipy.stats.hypergeom.sf(k - 1, M, n, N)),
    multiple_testing_correction: Callable[[dict[Hashable, float]], dict[
        Hashable, float]] = correction.benjamini_yekutieli,
) -> dict[nx.Graph, dict[int, dict[str, tuple[float, frozenset[str]]]]]:
    """
    Test communities for enrichment of large protein-specific measurements for
    each time of measurement and type of post-translational modification with
    respect to the protein-protein interaction network.

    Args:
        network: The protein-protein interaction network.
        communities: The communities of the protein-protein interaction network.
        measurement_ranges: Modification-specific range of averaged measurements
            by which proteins are categorized.
        measurement_score: Modification-specific functions used convert the
            measurement range bounds to their binary logarithm.
        site_average: Optional modification-specific functions to derive
            protein-specific measurements from site-specific measurements.
        replicate_average: Optional modification-specific functions to derive
            site-specific measurements from replicate measurements.
        enrichment_test: The statistical test used to assess enrichment of
            proteins with large measurements by a community.
        multiple_testing_correction: The procedure to correct for multiple
            testing of multiple communities, times of measurement and types of
            post-translational modification.

    Returns:
        Corrected p-values for the enrichment of protein-specific measurements
        by each community for each time of measurement and type of
        post-translational modification and associated proteins. If measurements
        are not averaged, proteins with any measurement exceeding the range are
        returned.
    """
    # Test the enrichment of average measurements a at most a lower or at least
    # an upper threshold by communities of the protein-protein interaction
    # network for different types of post-translational modification at
    # different times of measurement.
    p_values: dict[Hashable, float] = {}

    proteins = {}
    for time in get_times(network):
        for modification in get_modifications(network, time):
            # Determine the lower and upper threshold of the average measurement
            # for a type of post-translational modification at a time of
            # measurement.
            measurement_range = (
                measurement_score.get(
                    modification, lambda measurement, _: measurement)(
                        measurement_ranges.get(modification, (-1.0, 1.0))[0],
                        get_measurements(network, time, modification,
                                         site_average.get(modification),
                                         replicate_average.get(modification))),
                measurement_score.get(
                    modification, lambda measurement, _: measurement)(
                        measurement_ranges.get(modification, (-1.0, 1.0))[1],
                        get_measurements(network, time, modification,
                                         site_average.get(modification),
                                         replicate_average.get(modification))))

            # Determine the number of average measurements in the
            # protein-protein interaction network associated with measurements
            # for a type of post-translational modification at a time of
            # measurement.
            measurements = len([
                measurement for measurement in get_measurements(
                    network, time, modification, site_average.get(modification),
                    replicate_average.get(modification)) if measurement
            ])

            # Compile a map from communities of the protein-protein interaction
            # network to the number of average measurements from that community
            # associated with measurements for a type of post-translational
            # modification at a time of measurement.
            community_measurements = {
                community: len([
                    measurement for measurement in get_measurements(
                        community, time, modification,
                        site_average.get(modification),
                        replicate_average.get(modification)) if measurement
                ]) for community in communities
            }

            # Determine the number of average measurements in the
            # protein-protein interaction network associated with measurements
            # for a type of post-translational modification at a time of
            # measurement at most the lower or at least the upper threshold.
            target_measurements = len([
                measurement for measurement in get_measurements(
                    network, time, modification, site_average.get(modification),
                    replicate_average.get(modification))
                if measurement <= measurement_range[0] or
                measurement >= measurement_range[1]
            ])

            # Compile a map from communities of the protein-protein interaction
            # network to the number of average measurements from that community
            # associated with measurements for a type of post-translational
            # modification at a time of measurement at most the lower or at
            # least the upper threshold.
            target_community_measurements = {
                community: len([
                    measurement for measurement in get_measurements(
                        community, time, modification,
                        site_average.get(modification),
                        replicate_average.get(modification))
                    if measurement <= measurement_range[0] or
                    measurement >= measurement_range[1]
                ]) for community in communities
            }

            # Compute p-values for enrichment of average measurements at most
            # the lower or at least the upper threshold for a type of
            # post-translational modification at a time of measurement.
            p_values.update({
                (community, time, modification):
                enrichment_test(target_community_measurements[community],
                                measurements, target_measurements,
                                community_measurements[community])
                for community in communities
            })

            # Determine proteins associated with average measurements at most
            # the lower or at least the upper threshold for a type of
            # post-translational modification at a time of measurement.
            proteins.update({(community, time, modification):
                             get_proteins(community, time, modification,
                                          site_average.get(modification),
                                          replicate_average.get(modification),
                                          measurement_range)
                             for community in communities})

    # Adjust p-values to correct for testing multiple communities of the
    # protein-protein interaction network, times of measurement and types of
    # post-translational modification.
    p_values = multiple_testing_correction(p_values)

    # Return the p-values and decisive proteins associated with average
    # measurements at most the lower or at least the upper threshold for
    # different types of post-translational modification and times of
    # measurement.
    return {
        community: {
            time: {
                modification: (p_values[(community, time, modification)],
                               proteins[(community, time, modification)])
                for modification in get_modifications(network, time)
            } for time in get_times(network)
        } for community in communities
    }


def get_location(
    network: nx.Graph,
    communities: Iterable[nx.Graph],
    site_average: dict[str, Optional[Callable[[Iterable[float]], float]]],
    replicate_average: dict[str, Optional[Callable[[Iterable[float]], float]]],
    location_test: Callable[
        [Collection[float], Collection[float]], float] = lambda x, y: float(
            scipy.stats.mannwhitneyu(
                x, y, use_continuity=False, method="asymptotic").pvalue),
    multiple_testing_correction: Callable[[dict[Hashable, float]], dict[
        Hashable, float]] = correction.benjamini_yekutieli,
) -> dict[nx.Graph, dict[int, dict[str, float]]]:
    """
    Test communities for equal location of average measurements relative to the
    remaining protein-protein interaction network for each time of measurement
    and type of post-translational modification.

    Args:
        network: The protein-protein interaction network.
        communities: The communities of the protein-protein interaction network.
        site_average: Optional modification-specific functions to derive
            protein-specific measurements from site-specific measurements.
        replicate_average: Optional modification-specific functions to derive
            site-specific measurements from replicate measurements.
        location_test: The statistical test used to assess location of
            proteins with large measurements by a community.
        multiple_testing_correction: The procedure to correct for multiple
            testing of multiple communities, times of measurement and types of
            post-translational modification.

    Returns:
        Corrected p-values for the equality of location of average measurements
        of communities relative to the remaining protein-protein interaction
        network for each time of measurement and type of post-translational
        modification.
    """
    # Test the equality of location of the average measurements of communities
    # relative to the remaining protein-protein interaction network for
    # different types of post-translational modification at different times of
    # measurement.
    p_values: dict[Hashable, float] = {}

    for time in get_times(network):
        for modification in get_modifications(network, time):
            # Compile the distributions of average measurements for a type of
            # post-translational modification at a time of measurement in each
            # community of the protein-protein interaction network.
            community_measurements = [
                get_measurements(community, time, modification,
                                 site_average.get(modification),
                                 replicate_average.get(modification))
                for community in communities
            ]

            # Compile the distributions of average measurements for a type of
            # post-translational modification at a time of measurement in the
            # remaining protein-protein interaction network relative to each
            # community.
            network_measurements = [
                get_measurements(
                    nx.union_all([c
                                  for c in communities
                                  if c != community]), time, modification,
                    site_average.get(modification),
                    replicate_average.get(modification))
                for community in communities
            ]

            # Compute p-values for the equality of location of average
            # measurements for a type of post-translational modification at a
            # time of measurement of communities relative to the remaining
            # protein-protein interaction network.
            p_values.update({(community, time, modification):
                             location_test(community_measurements[i],
                                           network_measurements[i])
                             for i, community in enumerate(communities)
                             if community_measurements[i]})

    # Adjust p-values to correct for testing multiple communities of the
    # protein-protein interaction network, times of measurement and types of
    # post-translational modification.
    p_values = multiple_testing_correction(p_values)

    # Return the p-values.
    return {
        community: {
            time: {
                modification: p_values[(community, time, modification)]
                for modification in get_modifications(network, time)
                if (community, time, modification) in p_values
            } for time in get_times(network)
        } for community in communities
    }


def export(network: nx.Graph, basename: str) -> Optional[str]:
    """
    Exports the protein-protein interaction network if possible without
    overwriting an existing file.

    Args:
        network: The protein-protein interaction network.
        basename: The base file name.

    Returns:
        The file the protein-protein interaction network was exported to if
        there is no naming conflict.
    """
    # Avoid overwriting an existing file.
    if os.path.isfile(f"{basename}.graphml"):
        return None

    # Export the element tree of the protein-protein interaction network.
    nx.write_graphml_xml(network,
                         f"{basename}.graphml",
                         named_key_ids=True,
                         infer_numeric_types=True)

    # Return the file name of the protein-protein interaction network.
    return f"{basename}.graphml"
