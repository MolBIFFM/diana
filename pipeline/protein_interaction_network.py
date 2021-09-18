import bisect
import concurrent.futures
import itertools
import math
import os
import statistics

import networkx as nx
import pandas as pd
import scipy.stats

from pipeline.modularization import louvain
from pipeline.utilities import correction, uniprot


class ProteinInteractionNetwork(nx.Graph):
    def __init__(self):
        super().__init__()

    def annotate_proteins(self):
        for accessions, gene_name, protein_name, _ in uniprot.get_swissprot_entries(
        ):
            for protein in tuple(self):
                if protein.split("-")[0] in accessions:
                    if "-" in protein and not protein.split(
                            "-")[1].isnumeric():
                        nx.relabel_nodes(self,
                                         {protein: protein.split("-")[0]},
                                         copy=False)
                        protein = protein.split("-")[0]

                    if accessions.index(protein.split("-")[0]) == 0:
                        self.nodes[protein]["gene"] = gene_name
                        self.nodes[protein]["protein"] = protein_name
                    else:
                        if "-" in protein:
                            primary_protein = "{}-{}".format(
                                accessions[0],
                                protein.split("-")[1])
                        else:
                            primary_protein = accessions[0]

                        nx.relabel_nodes(self, {protein: primary_protein},
                                         copy=False)
                        self.nodes[primary_protein]["gene"] = gene_name
                        self.nodes[primary_protein]["protein"] = protein_name

    def remove_unannotated_proteins(self):
        self.remove_nodes_from([
            node for node in self if not (self.nodes[node].get("gene")
                                          and self.nodes[node].get("protein"))
        ])

    def add_genes_from(self, genes, taxon_identifier=9606):
        for accessions, gene_name, protein_name, taxon in uniprot.get_swissprot_entries(
        ):
            if taxon == taxon_identifier and gene_name in genes:
                self.add_node(accessions[0])
                self.nodes[accessions[0]]["gene"] = gene_name
                self.nodes[accessions[0]]["protein"] = protein_name

    def add_genes_from_table(
        self,
        file_name,
        gene_accession_column,
        gene_accession_format,
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
                gene_accession_format(row[gene_accession_column])
            ])

        self.add_genes_from(genes, taxon_identifier)

    def add_proteins_from(self, proteins):
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
        for accessions, gene_name, protein_name, _ in uniprot.get_swissprot_entries(
        ):
            for i, accession in enumerate(accessions):
                if accession in proteins_isoform:
                    for isoform in proteins_isoform[accession]:
                        if isoform == "0":
                            self.add_node(accessions[0])
                            self.nodes[accessions[0]]["gene"] = gene_name
                            self.nodes[accessions[0]]["protein"] = protein_name
                        else:
                            self.add_node("{}-{}".format(
                                accessions[0], isoform))
                            self.nodes["{}-{}".format(
                                accessions[0], isoform)]["gene"] = gene_name
                            self.nodes["{}-{}".format(
                                accessions[0],
                                isoform)]["protein"] = protein_name

                    if i > 0:
                        if accession not in primary_accession:
                            primary_accession[accession] = set()

                        primary_accession[accession].add(accessions[0])

        return primary_accession

    def add_proteins_from_table(
        self,
        file_name,
        protein_accession_column,
        protein_accession_format=lambda entry: [entry],
        time=0,
        ptm="",
        position_column="",
        position_format=lambda entry: [entry],
        replicates=[],
        sheet_name=0,
        header=0,
        num_sites=1000,
        num_replicates=1,
        combine_replicates=statistics.mean,
        convert_measurement=math.log2,
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
                usecols=[protein_accession_column] +
                [column
                 for column in (position_column, ) if column] + replicates,
                dtype={
                    protein_accession_column: str,
                    **{
                        column: str
                        for column in (position_column, ) if column
                    },
                    **{column: float
                       for column in replicates},
                },
            )
        else:
            table = pd.read_table(
                file_name,
                header=header,
                usecols=[protein_accession_column] +
                [column
                 for column in (position_column, ) if column] + replicates,
                dtype={
                    protein_accession_column: str,
                    **{
                        column: str
                        for column in (position_column, ) if column
                    },
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
                protein_accession_format(row[protein_accession_column])
            ]

            if replicates:
                if position_column and not pd.isna(row[position_column]):
                    positions = [
                        int(position)
                        for position in position_format(row[position_column])
                    ]
                else:
                    positions = []

                if len(protein_accessions) != len(positions):
                    if len(protein_accessions) > len(positions):
                        positions.extend([
                            0 for _ in range(len(positions),
                                             len(protein_accessions))
                        ])
                    else:
                        positions = positions[:len(protein_accessions)]

                measurements = [
                    row[replicate] for replicate in replicates
                    if not pd.isna(row[replicate])
                ]

                if len(measurements) >= min(num_replicates, len(replicates)):
                    for protein_accession, position in zip(
                            protein_accessions, positions):
                        if protein_accession not in proteins:
                            proteins[protein_accession] = []

                        bisect.insort(
                            proteins[protein_accession],
                            (
                                position,
                                convert_measurement(
                                    combine_replicates(measurements)),
                            ),
                        )
            else:
                for protein_accession in protein_accessions:
                    if protein_accession not in proteins:
                        proteins[protein_accession] = []

        primary_accession = self.add_proteins_from(proteins.keys())

        for protein in tuple(proteins):
            if protein.split("-")[0] in primary_accession:
                for accession in sorted(
                        primary_accession[protein.split("-")[0]]):
                    if "-" in protein:
                        primary_protein = "{}-{}".format(
                            accession,
                            protein.split("-")[1])
                    else:
                        primary_protein = accession

                    if primary_protein in proteins:
                        for mi in range(len(proteins[protein])):
                            for mj in range(len(proteins[primary_protein])):
                                if (proteins[protein][mi][1] ==
                                        proteins[primary_protein][mj][1]):
                                    proteins[protein][mi] = tuple()
                                    break

                        for measurement in proteins[protein]:
                            if measurement:
                                bisect.insort(
                                    proteins[primary_protein],
                                    measurement,
                                )

                    else:
                        proteins[primary_protein] = proteins[protein]

                del proteins[protein]

        for protein in self:
            if protein not in proteins:
                continue

            proteins[protein] = sorted(
                sorted(
                    proteins[protein],
                    key=lambda item: abs(item[1]),
                )[-num_sites:])

            for i in range(len(proteins[protein])):
                self.nodes[protein]["{} {} {}".format(
                    time, ptm, i + 1)] = proteins[protein][i][1]

    def get_times(self):
        return tuple(
            sorted(
                set(
                    int(change.split(" ")[0]) for protein in self
                    for change in self.nodes[protein]
                    if len(change.split(" ")) == 3
                    and change.split(" ")[0].isnumeric())))

    def get_post_translational_modifications(self, time):
        return tuple(
            sorted(
                set(
                    change.split(" ")[1] for protein in self
                    for change in self.nodes[protein]
                    if len(change.split(" ")) == 3
                    and change.split(" ")[0] == str(time))))

    def get_sites(self, time, ptm):
        return max(
            int(change.split(" ")[2]) for protein in self
            for change in self.nodes[protein]
            if len(change.split(" ")) == 3 and change.split(" ")[0] == str(
                time) and change.split(" ")[1] == ptm)

    def set_post_translational_modification(self):
        for time in self.get_times():
            for protein in self:
                self.nodes[protein][
                    "post-translational modification {}".format(
                        time)] = " ".join(
                            sorted(
                                set(
                                    change.split(" ")[1]
                                    for change in self.nodes[protein]
                                    if len(change.split(" ")) == 3
                                    and change.split(" ")[0] == str(time))))

    def get_changes(self,
                    time,
                    ptm,
                    combine_sites=lambda sites: max(sites, key=abs)):
        changes = []
        for protein in self:
            sites = [
                self.nodes[protein][change] for change in self.nodes[protein]
                if len(change.split(" ")) == 3 and change.split(" ")[0] == str(
                    time) and change.split(" ")[1] == ptm
            ]

            if sites:
                changes.append(combine_sites(sites))

        return changes

    def get_z_score_range(
            self,
            time,
            modification,
            thresholds,
            combine_sites=lambda sites: max(sites, key=abs),
    ):
        changes = self.get_changes(time, modification, combine_sites)
        mean = statistics.mean(changes)
        stdev = statistics.stdev(changes, xbar=mean)
        return (thresholds[0] * stdev + mean, thresholds[1] * stdev + mean)

    def get_z_score(
            self,
            time,
            modification,
            change,
            combine_sites=lambda sites: max(sites, key=abs),
    ):
        changes = self.get_changes(time, modification, combine_sites)
        mean = statistics.mean(changes)
        stdev = statistics.stdev(changes, xbar=mean)
        return (change - mean) / stdev

    def get_proportion_range(
            self,
            time,
            modification,
            thresholds,
            combine_sites=lambda sites: max(sites, key=abs),
    ):
        changes = sorted(self.get_changes(time, modification, combine_sites))
        proportion_range = [0.0, 0.0]

        for i in range(len(changes)):
            if changes[i] == changes[i + 1]:
                continue

            if i / len(changes) > thresholds[0]:
                proportion_range[0] = changes[i - 1]
                break

        for i in range(len(changes) - 1, -1, -1):
            if changes[i] == changes[i - 1]:
                continue

            if (len(changes) - i) / len(changes) > 1.0 - thresholds[1]:
                proportion_range[1] = changes[i + 1]
                break

        return tuple(proportion_range)

    def get_proportion(
            self,
            time,
            modification,
            change,
            combine_sites=lambda sites: max(sites, key=abs),
    ):
        changes = self.get_changes(time, modification, combine_sites)
        return (min(
            len([c for c in changes if c <= change]),
            len([c for c in changes if c >= change]),
        ) / len(changes))

    def set_changes(
        self,
        combine_sites=lambda sites: max(sites, key=abs),
        changes=(-1.0, 1.0),
        get_range=lambda time, modification, changes, combine_sites: changes,
    ):
        times = self.get_times()
        modifications = {
            time: self.get_post_translational_modifications(time)
            for time in times
        }

        mid_range = {
            time: {
                modification: get_range(
                    time,
                    modification,
                    changes,
                    combine_sites,
                )
                for modification in modifications[time]
            }
            for time in times
        }

        for time in times:
            for protein in self:
                ptm = {}
                for modification in modifications[time]:
                    sites = [
                        self.nodes[protein][change]
                        for change in self.nodes[protein]
                        if len(change.split(" ")) == 3 and change.split(" ")[0]
                        == str(time) and change.split(" ")[1] == modification
                    ]

                    if sites:
                        combined_change = combine_sites(sites)

                        if combined_change >= 1.0 * mid_range[time][
                                modification][1]:
                            ptm[modification] = "up"
                        elif 1.0 * mid_range[time][modification][
                                1] > combined_change >= 0.5 * mid_range[time][
                                    modification][1]:
                            ptm[modification] = "mid_up"
                        elif 0.5 * mid_range[time][modification][
                                1] > combined_change > 0.5 * mid_range[time][
                                    modification][0]:
                            ptm[modification] = "mid"
                        elif 0.5 * mid_range[time][modification][
                                0] >= combined_change > 1.0 * mid_range[time][
                                    modification][0]:
                            ptm[modification] = "mid_down"
                        else:
                            ptm[modification] = "down"

                if ptm:
                    if set(ptm.values()) == {"up"}:
                        self.nodes[protein]["change {}".format(time)] = "up"
                    elif set(ptm.values()) == {"mid_up", "up"}:
                        self.nodes[protein]["change {}".format(time)] = "up"
                    elif set(ptm.values()) == {"mid", "mid_up", "up"}:
                        self.nodes[protein]["change {}".format(time)] = "up"

                    elif set(ptm.values()) == {"mid_up"}:
                        self.nodes[protein]["change {}".format(
                            time)] = "mid_up"
                    elif set(ptm.values()) == {"mid", "mid_up"}:
                        self.nodes[protein]["change {}".format(
                            time)] = "mid_up"

                    elif set(ptm.values()) == {"mid"}:
                        self.nodes[protein]["change {}".format(time)] = "mid"

                    elif set(ptm.values()) == {"mid", "mid_down"}:
                        self.nodes[protein]["change {}".format(
                            time)] = "mid_down"
                    elif set(ptm.values()) == {"mid_down"}:
                        self.nodes[protein]["change {}".format(
                            time)] = "mid_down"

                    elif set(ptm.values()) == {"mid", "mid_down", "down"}:
                        self.nodes[protein]["change {}".format(time)] = "down"
                    elif set(ptm.values()) == {"mid_down", "down"}:
                        self.nodes[protein]["change {}".format(time)] = "down"
                    elif set(ptm.values()) == {"down"}:
                        self.nodes[protein]["change {}".format(time)] = "down"
                    else:
                        self.nodes[protein]["change {}".format(
                            time)] = " ".join(
                                "{}_{}".format(modification, ptm[modification])
                                for modification in sorted(ptm.keys()))

                else:
                    self.nodes[protein]["change {}".format(time)] = "mid"

    def get_databases(self):
        return tuple(
            sorted(
                set(database for edge in self.edges
                    for database in self.edges[edge]).intersection(
                        {"BioGRID", "CORUM", "IntAct", "Reactome", "STRING"})))

    def set_edge_weights(
        self,
        weight=lambda confidence_scores: int(bool(confidence_scores)),
        attribute="weight",
    ):
        databases = self.get_databases()

        for edge in self.edges:
            self.edges[edge][attribute] = weight({
                database: self.edges[edge][database]
                for database in databases if database in self.edges[edge]
            })

    def get_modules(
        self,
        module_size=0,
        combine_module_sizes=statistics.mean,
        weight="weight",
        algorithm=louvain.louvain,
    ):
        G = self.copy()
        G.remove_nodes_from(tuple(nx.isolates(G)))

        if G.number_of_nodes() == 0:
            return []

        communities = algorithm(G, weight=weight)

        if not module_size:
            return communities

        while (combine_module_sizes(
            [len(community) for community in communities]) > module_size):
            max_community_size = max(
                len(community) for community in communities)

            indices = [
                communities.index(community) for community in communities
                if len(community) == max_community_size
            ]

            subdivision = False

            with concurrent.futures.ProcessPoolExecutor() as executor:
                for i, subdivided_community in zip(
                        indices,
                        executor.map(
                            algorithm,
                            (G.subgraph(communities[j]) for j in indices),
                            itertools.repeat(weight),
                        ),
                ):
                    if len(subdivided_community) > 1:
                        communities[i:i + 1] = subdivided_community

                        indices[indices.index(i) + 1:] = [
                            k + len(subdivided_community) - 1
                            for k in indices[indices.index(i) + 1:]
                        ]

                        subdivision = True

            if not subdivision:
                break

        return [community for community in communities if len(community) > 1]

    def get_proteins(
            self,
            time,
            ptm,
            change_filter=lambda combined_sites: bool(combined_sites),
            combine_sites=lambda sites: max(sites, key=abs),
    ):
        proteins = []
        for protein in self:
            sites = [
                self.nodes[protein][change] for change in self.nodes[protein]
                if len(change.split(" ")) == 3 and change.split(" ")[0] == str(
                    time) and change.split(" ")[1] == ptm
            ]

            if sites and change_filter(combine_sites(sites)):
                proteins.append(protein)

        return proteins

    def get_module_change_enrichment(
        self,
        p=0.05,
        changes=(-1.0, 1.0),
        get_range=lambda time, ptm, changes, combine_sites: changes,
        combine_sites=lambda sites: max(sites, key=abs),
        module_size=0,
        combine_module_sizes=statistics.mean,
        weight="weight",
        test="two-sided",
        algorithm=louvain.louvain,
    ):
        modules = {
            i: module
            for i, module in enumerate(
                self.get_modules(module_size,
                                 combine_module_sizes,
                                 weight=weight,
                                 algorithm=algorithm))
        }

        p_values = {}
        p_adjusted = {}
        modules_filtered = {}

        for time in self.get_times():
            for ptm in self.get_post_translational_modifications(time):
                M = len(self.get_proteins(time, ptm))
                N = [
                    len(self.subgraph(modules[i]).get_proteins(time, ptm))
                    for i in modules
                ]

                mid_range = get_range(
                    time,
                    ptm,
                    changes,
                    combine_sites,
                )

                if test == "one-sided positive":
                    n = len(
                        self.get_proteins(
                            time,
                            ptm,
                            lambda change: change > mid_range[1],
                            combine_sites=combine_sites,
                        ))

                    k = [
                        len(
                            self.subgraph(modules[i]).get_proteins(
                                time,
                                ptm,
                                lambda change: change > mid_range[1],
                                combine_sites=combine_sites,
                            )) for i in modules
                    ]

                elif test == "one-sided negative":
                    n = len(
                        self.get_proteins(
                            time,
                            ptm,
                            lambda change: change < mid_range[0],
                            combine_sites=combine_sites,
                        ))

                    k = [
                        len(
                            self.subgraph(modules[i]).get_proteins(
                                time,
                                ptm,
                                lambda change: change < mid_range[0],
                                combine_sites=combine_sites,
                            )) for i in modules
                    ]

                elif test == "two-sided":
                    n = len(
                        self.get_proteins(
                            time,
                            ptm,
                            lambda change: change < mid_range[0],
                            combine_sites=combine_sites,
                        )) + len(
                            self.get_proteins(
                                time,
                                ptm,
                                lambda change: change > mid_range[1],
                                combine_sites=combine_sites,
                            ))

                    k = [
                        len(
                            self.subgraph(modules[i]).get_proteins(
                                time,
                                ptm,
                                lambda change: change < mid_range[0],
                                combine_sites=combine_sites,
                            )) + len(
                                self.subgraph(modules[i]).get_proteins(
                                    time,
                                    ptm,
                                    lambda change: change > mid_range[1],
                                    combine_sites=combine_sites,
                                )) for i in modules
                    ]

                else:
                    n = 0
                    k = [0 for _ in modules]

                p_values = {
                    i: scipy.stats.hypergeom.sf(k[i] - 1, M, n, N[i])
                    for i in range(len(modules))
                }

                for module, p_value in correction.benjamini_hochberg(
                        p_values).items():
                    if p_value <= p:
                        if time not in p_adjusted:
                            p_adjusted[time] = {}
                        if ptm not in p_adjusted[time]:
                            p_adjusted[time][ptm] = {}

                        modules_filtered[module] = modules[module]
                        p_adjusted[time][ptm][module] = p_value

        return modules_filtered, p_adjusted

    def get_neighborhood(self, protein, k=1, isoforms=True):
        if isoforms:
            nodes = {
                node
                for node in self if node.split("-")[0] == protein.split("-")[0]
            }
        else:
            nodes = {protein.split("-")[0]}

        for _ in range(k):
            nodes.update(
                set.union(*(set(self.neighbors(node)) for node in nodes)))

        return self.subgraph(nodes)
