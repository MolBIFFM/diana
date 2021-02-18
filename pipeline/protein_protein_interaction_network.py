import bisect
import itertools
import json
import math
import statistics

import networkx as nx
import scipy.stats
import pandas as pd

from pipeline.configuration import data
from pipeline.utilities import download, mitab, modularization


class ProteinProteinInteractionNetwork(nx.Graph):
    def __init__(self):
        super().__init__()

    def add_proteins_from_spreadsheet(
        self,
        file_name,
        ptm,
        time,
        protein_accession_column,
        position_col,
        replicates,
        protein_accession_format,
        position_format,
        sheet_name=0,
        header=0,
        num_sites=1,
        num_replicates=1,
        merge_replicates=statistics.mean,
        convert_measurement=math.log2,
    ):
        proteins = {}
        for _, row in pd.read_excel(
            file_name,
            sheet_name=sheet_name,
            header=header,
            usecols=[protein_accession_column, position_col] + replicates,
            dtype={
                protein_accession_column: str,
                position_col: str,
                **{replicate: float for replicate in replicates},
            },
        ).iterrows():

            if pd.isna(row[protein_accession_column]) or pd.isna(row[position_col]):
                continue

            measurements = [row[repl] for repl in replicates if not pd.isna(row[repl])]

            if len(measurements) >= min(num_replicates, len(replicates)):
                protein_accessions = [
                    str(protein_accession)
                    for protein_accession in protein_accession_format(
                        row[protein_accession_column]
                    )
                ]
                positions = [
                    int(position) for position in position_format(row[position_col])
                ]

                if len(protein_accessions) >= len(positions):
                    positions.extend([0 for _ in range(len(positions), len(proteins))])

                    for protein_accession, position in zip(
                        protein_accessions, positions
                    ):
                        if "-" in protein_accession:
                            protein, isoform = protein_accession.split("-")
                        else:
                            protein, isoform = protein_accession, "1"

                        if protein not in proteins:
                            proteins[protein] = {}
                        if isoform not in proteins[protein]:
                            proteins[protein][isoform] = []

                        bisect.insort(
                            proteins[protein][isoform],
                            (
                                position,
                                convert_measurement(merge_replicates(measurements)),
                            ),
                        )
                else:
                    for protein in proteins:
                        for position in position:
                            if "-" in protein_accession:
                                protein, isoform = protein_accession.split("-")
                            else:
                                protein, isoform = protein_accession, "1"

                            if protein not in proteins:
                                proteins[protein] = {}
                            if isoform not in proteins[protein]:
                                proteins[protein][isoform] = []

                            bisect.insort(
                                proteins[protein][isoform],
                                (
                                    position,
                                    convert_measurement(merge_replicates(measurements)),
                                ),
                            )

        reviewed_proteins, primary_accession, gene_name, protein_name = {}, {}, {}, {}
        accessions, gene_names, protein_names = [], {}, {}
        recommended_name = False
        for line in download.iterate_txt(data.UNIPROT_SWISS_PROT):
            if line.split("   ")[0] == "AC":
                accessions.extend(line.split("   ")[1].rstrip(";").split("; "))

            elif line.split("   ")[0] == "GN":
                for entry in line.split("   ")[1].rstrip(";").split("; "):
                    if "=" in entry:
                        gene_names[entry.split("=")[0]] = (
                            entry.split("=")[1].split("{")[0].rstrip()
                        )

            elif line.split("   ")[0] == "DE":
                if line.split(" ", 1)[1].lstrip().split(":", 1)[0] == "RecName":
                    entries = line.split(":", 1)[1].lstrip().rstrip(";").split("; ")
                    recommended_name = True
                elif line.split(" ", 1)[1].lstrip().split(":", 1)[0] == "AltName":
                    entries = []
                    recommended_name = False
                elif line.split(" ", 1)[1].lstrip().split(":", 1)[0] in (
                    "Flags",
                    "Contains",
                ):
                    entries = []
                elif recommended_name:
                    entries = line.split(" ", 1)[1].lstrip().rstrip(";").split("; ")

                for entry in entries:
                    if "=" in entry:
                        protein_names[entry.split("=")[0]] = (
                            entry.split("=")[1].split("{")[0].rstrip()
                        )

            elif line == "//":
                for i, accession in enumerate(accessions):
                    if accession in proteins:
                        reviewed_proteins[accession] = proteins[accession]
                        gene_name[accession] = gene_names.get("Name", "NA")
                        protein_name[accession] = protein_names.get("Full", "NA")
                        if i > 0:
                            primary_accession[accession] = accessions[0]
                accessions.clear()
                gene_names.clear()
                protein_names.clear()

        proteins = reviewed_proteins

        for protein in primary_accession:
            if primary_accession[protein] in proteins:
                for primary_isoform in list(proteins[primary_accession[protein]]):
                    primary_positions = {}
                    for position, measurement in proteins[primary_accession[protein]][
                        primary_isoform
                    ]:
                        if position not in primary_positions:
                            primary_positions[position] = set()
                        primary_positions[position].add(measurement)

                    for secondary_isoform in list(proteins[protein]):
                        secondary_positions = {}
                        for position, measurement in proteins[protein][
                            secondary_isoform
                        ]:
                            if position not in secondary_positions:
                                secondary_positions[position] = set()
                            secondary_positions[position].add(measurement)

                        if [
                            measurements
                            for position, measurements in sorted(
                                secondary_positions.items()
                            )
                        ] == [
                            measurements
                            for position, measurements in sorted(
                                primary_positions.items()
                            )
                        ]:
                            positions = [
                                secondary_position
                                if not primary_position
                                else primary_position
                                for primary_position, secondary_position in zip(
                                    sorted(primary_positions.keys()),
                                    sorted(secondary_positions.keys()),
                                )
                            ]

                        proteins[primary_accession[protein]][primary_isoform] = sorted(
                            zip(
                                positions,
                                [
                                    measurement
                                    for position, measurement in proteins[
                                        primary_accession[protein]
                                    ][primary_isoform]
                                ],
                            )
                        )

                        del proteins[protein][secondary_isoform]
                        break

                if not proteins[protein]:
                    del proteins[protein]
                    if protein in gene_name:
                        del gene_name[protein]
            else:
                proteins[primary_accession[protein]] = proteins.pop(protein)
                gene_name[primary_accession[protein]] = gene_name.pop(protein)
                protein_name[primary_accession[protein]] = protein_name.pop(protein)

        for protein in proteins:
            for isoform in proteins[protein]:
                if isoform != "1":
                    self.add_node("-".join([protein, isoform]))
                    self.nodes["-".join([protein, isoform])]["gene name"] = gene_name[
                        protein
                    ]
                    self.nodes["-".join([protein, isoform])][
                        "protein name"
                    ] = protein_name[protein]
                else:
                    self.add_node(protein)
                    self.nodes[protein]["gene name"] = gene_name[protein]
                    self.nodes[protein]["protein name"] = protein_name[protein]

                if len(proteins[protein][isoform]) > num_sites:
                    proteins[protein][isoform] = sorted(
                        sorted(
                            proteins[protein][isoform],
                            key=lambda tp: tp[1],
                            reverse=True,
                        )[:num_sites]
                    )
                for i in range(len(proteins[protein][isoform])):
                    if isoform != "1":
                        self.nodes["-".join([protein, isoform])][
                            "{} {} {}".format(time, ptm, i + 1)
                        ] = proteins[protein][isoform][i][1]
                    else:
                        self.nodes[protein][
                            "{} {} {}".format(time, ptm, i + 1)
                        ] = proteins[protein][isoform][i][1]

                yield [
                    gene_name[protein],
                    "-".join([protein, isoform]) if isoform != "1" else protein,
                    protein_name[protein],
                    [
                        proteins[protein][isoform][i][1]
                        for i in range(len(proteins[protein][isoform]))
                    ],
                ]

    def get_times(self):
        return sorted(
            set(
                int(change.split(" ")[0])
                for protein in self
                for change in self.nodes[protein]
                if len(change.split(" ")) == 3 and change.split(" ")[0].isnumeric()
            )
        )

    def get_post_translational_modifications(self, time):
        return sorted(
            set(
                change.split(" ")[1]
                for protein in self
                for change in self.nodes[protein]
                if len(change.split(" ")) == 3
                and change.split(" ")[0].isnumeric()
                and change.split(" ")[0] == str(time)
            )
        )

    def get_sites(self, time, ptm):
        return max(
            int(change.split(" ")[2])
            for protein in self
            for change in self.nodes[protein]
            if len(change.split(" ")) == 3
            and change.split(" ")[0].isnumeric()
            and change.split(" ")[0] == str(time)
            and change.split(" ")[1] == ptm
        )

    def set_post_translational_modification_data(self):
        for time in self.get_times():
            for protein in self:
                self.nodes[protein][
                    "post-translational modification {}".format(time)
                ] = " ".join(
                    sorted(
                        set(
                            change.split(" ")[1]
                            for change in self.nodes[protein]
                            if len(change.split(" ")) == 3
                            and change.split(" ")[0] == str(time)
                        )
                    )
                )

    def get_changes(self, time, ptm, merge_sites=lambda sites: max(sites, key=abs)):
        changes = []
        for protein in self:
            sites = [
                self.nodes[protein][change]
                for change in self.nodes[protein]
                if len(change.split(" ")) == 3
                and change.split(" ")[0].isnumeric()
                and change.split(" ")[0] == str(time)
                and change.split(" ")[1] == ptm
            ]

            if sites:
                changes.append(merge_sites(sites))

        return changes

    def get_range_from_z_score(
        self,
        time,
        ptm,
        merge_sites=lambda sites: max(sites, key=abs),
        threshold=(-2.0, 2.0),
    ):
        changes = self.get_changes(time, ptm, merge_sites)
        mean = statistics.mean(changes)
        std = statistics.pstdev(changes, mu=mean)

        return (threshold[0] * std + mean, threshold[1] * std + mean)

    def set_change_data(
        self,
        merge_sites=lambda sites: max(sites, key=abs),
        mid_range_thresholds=(-1.0, 1.0),
        get_mid_range=lambda time, ptm, merge_sites, mid_range_thresholds: mid_range_thresholds,
    ):
        times = self.get_times()
        post_translational_modifications = {
            time: self.get_post_translational_modifications(time) for time in times
        }

        mid_range = {
            time: {
                post_translational_modification: get_mid_range(
                    time,
                    post_translational_modification,
                    merge_sites,
                    mid_range_thresholds,
                )
                for post_translational_modification in post_translational_modifications[
                    time
                ]
            }
            for time in times
        }

        for time in times:
            for protein in self:
                ptms = {}
                for post_translational_modification in post_translational_modifications[
                    time
                ]:
                    sites = [
                        self.nodes[protein][change]
                        for change in self.nodes[protein]
                        if len(change.split(" ")) == 3
                        and change.split(" ")[0].isnumeric()
                        and change.split(" ")[0] == str(time)
                        and change.split(" ")[1] == post_translational_modification
                    ]
                    if sites:
                        ptms[post_translational_modification] = merge_sites(sites)

                if ptms:
                    if all(
                        change
                        >= 0.5 * mid_range[time][post_translational_modification][1]
                        for change in ptms.values()
                    ):
                        if any(
                            change
                            >= mid_range[time][post_translational_modification][1]
                            for change in ptms.values()
                        ):
                            self.nodes[protein]["change {}".format(time)] = "up"
                        else:
                            self.nodes[protein]["change {}".format(time)] = "mid up"

                    elif all(
                        change
                        <= 0.5 * mid_range[time][post_translational_modification][0]
                        for change in ptms.values()
                    ):
                        if any(
                            change
                            <= mid_range[time][post_translational_modification][0]
                            for change in ptms.values()
                        ):
                            self.nodes[protein]["change {}".format(time)] = "down"
                        else:
                            self.nodes[protein]["change {}".format(time)] = "mid down"

                    elif all(
                        change
                        <= 0.5 * mid_range[time][post_translational_modification][0]
                        or change
                        >= 0.5 * mid_range[time][post_translational_modification][1]
                        for change in ptms.values()
                    ):
                        self.nodes[protein]["change {}".format(time)] = " ".join(
                            sorted(
                                [
                                    "{} up".format(ptm)
                                    if change > 0.0
                                    else "{} down".format(ptm)
                                    for ptm, change in ptms.items()
                                ]
                            )
                        )
                    else:
                        self.nodes[protein]["change {}".format(time)] = "mid"

    def get_databases(self):
        return sorted(
            set(
                database for edge in self.edges for database in self.edges[edge]
            ).intersection({"BioGRID", "IntAct", "STRING", "Reactome"})
        )

    def set_edge_weights(
        self,
        weight=lambda confidence_scores: bool(confidence_scores),
        attribute="weight",
    ):
        databases = self.get_databases()

        for edge in self.edges:
            self.edges[edge][attribute] = weight(
                {
                    database: self.edges[edge][database]
                    for database in databases
                    if database in self.edges[edge]
                }
            )

    def add_interactions_from_biogrid(
        self, experimental_system=[], multi_validated_physical=False
    ):
        uniprot = {}
        for row in download.iterate_tabular_txt(
            data.UNIPROT_ID_MAP,
            delimiter="\t",
            usecols=[0, 1, 2],
        ):
            if row[1] == "BioGRID" and row[0] in self.nodes:
                uniprot[int(row[2])] = row[0]

        for row in download.iterate_tabular_txt(
            data.BIOGRID_ID_MAP_ZIP_ARCHIVE,
            zip_file=data.BIOGRID_ID_MAP,
            delimiter="\t",
            header=20,
            usecols=[
                "BIOGRID_ID",
                "IDENTIFIER_VALUE",
                "IDENTIFIER_TYPE",
                "ORGANISM_OFFICIAL_NAME",
            ],
        ):
            if (
                row["IDENTIFIER_TYPE"] in ("UNIPROT-ACCESSION", "UNIPROT-ISOFORM")
                and row["ORGANISM_OFFICIAL_NAME"] == "Homo sapiens"
                and row["IDENTIFIER_VALUE"] in self.nodes
            ):
                uniprot[int(row["BIOGRID_ID"])] = row["IDENTIFIER_VALUE"]

        for row in download.iterate_tabular_txt(
            data.BIOGRID_MV_PHYSICAL_ZIP_ARCHIVE
            if multi_validated_physical
            else data.BIOGRID_ZIP_ARCHIVE,
            zip_file=data.BIOGRID_MV_PHYSICAL
            if multi_validated_physical
            else data.BIOGRID,
            delimiter="\t",
            header=0,
            usecols=[
                "BioGRID ID Interactor A",
                "BioGRID ID Interactor B",
                "Experimental System",
                "Experimental System Type",
                "Organism ID Interactor A",
                "Organism ID Interactor B",
                "Throughput",
            ],
        ):
            if (
                uniprot.get(row["BioGRID ID Interactor A"])
                and uniprot.get(row["BioGRID ID Interactor B"])
                and uniprot[row["BioGRID ID Interactor A"]]
                != uniprot[row["BioGRID ID Interactor B"]]
                and (
                    not experimental_system
                    or row["Experimental System"] in experimental_system
                )
                and row["Experimental System Type"] == "physical"
                and row["Organism ID Interactor A"] == 9606
                and row["Organism ID Interactor B"] == 9606
            ):
                self.add_edge(
                    uniprot[row["BioGRID ID Interactor A"]],
                    uniprot[row["BioGRID ID Interactor B"]],
                )
                self.edges[
                    uniprot[row["BioGRID ID Interactor A"]],
                    uniprot[row["BioGRID ID Interactor B"]],
                ]["BioGRID"] = 0.5

                yield (
                    uniprot[row["BioGRID ID Interactor A"]],
                    uniprot[row["BioGRID ID Interactor B"]],
                    self.edges[
                        uniprot[row["BioGRID ID Interactor A"]],
                        uniprot[row["BioGRID ID Interactor B"]],
                    ]["BioGRID"],
                )

    def add_interactions_from_intact(
        self, interaction_detection_methods=[], interaction_types=[], mi_score=0.0
    ):
        for row in download.iterate_tabular_txt(
            data.INTACT_ZIP_ARCHIVE,
            zip_file=data.INTACT,
            delimiter="\t",
            header=0,
            usecols=[
                "#ID(s) interactor A",
                "ID(s) interactor B",
                "Alt. ID(s) interactor A",
                "Alt. ID(s) interactor B",
                "Taxid interactor A",
                "Taxid interactor B",
                "Interaction detection method(s)",
                "Interaction type(s)",
                "Confidence value(s)",
            ],
        ):

            if not (
                interactor_a := mitab.get_id_from_namespace(
                    row["#ID(s) interactor A"], "uniprotkb"
                )
            ):
                if not (
                    interactor_a := mitab.get_id_from_namespace(
                        row["Alt. ID(s) interactor A"], "uniprotkb"
                    )
                ):
                    continue
            if interactor_a not in self:
                continue

            if not (
                interactor_b := mitab.get_id_from_namespace(
                    row["ID(s) interactor B"], "uniprotkb"
                )
            ):
                if not (
                    interactor_b := mitab.get_id_from_namespace(
                        row["Alt. ID(s) interactor B"], "uniprotkb"
                    )
                ):
                    continue
            if interactor_b not in self:
                continue

            if interactor_a == interactor_b:
                continue

            if not (
                mitab.namespace_has_id(row["Taxid interactor A"], "taxid", "9606")
                and mitab.namespace_has_id(row["Taxid interactor B"], "taxid", "9606")
            ):
                continue

            if not (
                not interaction_detection_methods
                or mitab.namespace_has_any_term_from(
                    row["Interaction detection method(s)"],
                    "psi-mi",
                    interaction_detection_methods,
                )
            ):
                continue

            if not (
                not interaction_types
                or mitab.namespace_has_any_term_from(
                    row["Interaction type(s)"], "psi-mi", interaction_types
                )
            ):
                continue

            if score := mitab.get_id_from_namespace(
                row["Confidence value(s)"], "intact-miscore"
            ):
                score = float(score)
                if score < mi_score:
                    continue
            else:
                continue

            if self.has_edge(interactor_a, interactor_b):
                self.edges[interactor_a, interactor_b]["IntAct"] = max(
                    score, self.edges[interactor_a, interactor_b].get("IntAct", 0.0)
                )
            else:
                self.add_edge(interactor_a, interactor_b)
                self.edges[interactor_a, interactor_b]["IntAct"] = score

            yield (
                interactor_a,
                interactor_b,
                self.edges[interactor_a, interactor_b]["IntAct"],
            )

    def add_interactions_from_string(
        self,
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
    ):

        uniprot = {}
        for row in download.iterate_tabular_txt(
            data.UNIPROT_ID_MAP,
            delimiter="\t",
            usecols=[0, 1, 2],
        ):
            if row[1] == "STRING" and row[0] in self.nodes:
                uniprot[row[2]] = row[0]

        for row in download.iterate_tabular_txt(data.STRING_ID_MAP, usecols=[1, 2]):
            if row[1].split("|")[0] in self.nodes:
                uniprot[row[2]] = row[1].split("|")[0]

        thresholds = {
            column: threshold
            for column, threshold in {
                "neighborhood": neighborhood,
                "neighborhood_transferred": neighborhood_transferred,
                "fusion": fusion,
                "cooccurence": cooccurence,
                "homology": homology,
                "coexpression": coexpression,
                "coexpression_transferred": coexpression_transferred,
                "experiments": experiments,
                "experiments_transferred": experiments_transferred,
                "database": database,
                "database_transferred": database_transferred,
                "textmining": textmining,
                "textmining_transferred": textmining_transferred,
            }.items()
            if threshold
        }
        thresholds["combined_score"] = combined_score

        for row in download.iterate_tabular_txt(
            data.STRING_PHYSICAL if physical else data.STRING,
            delimiter=" ",
            header=0,
            usecols=["protein1", "protein2"] + list(thresholds.keys()),
        ):
            if (
                uniprot.get(row["protein1"])
                and uniprot.get(row["protein2"])
                and uniprot[row["protein1"]] != uniprot[row["protein2"]]
                and all(
                    row[column] / 1000 >= thresholds[column] for column in thresholds
                )
            ):
                if self.has_edge(uniprot[row["protein1"]], uniprot[row["protein2"]]):
                    self.edges[uniprot[row["protein1"]], uniprot[row["protein2"]]][
                        "STRING"
                    ] = max(
                        row["combined_score"] / 1000,
                        self.edges[
                            uniprot[row["protein1"]], uniprot[row["protein2"]]
                        ].get("STRING", 0.0),
                    )
                else:
                    self.add_edge(uniprot[row["protein1"]], uniprot[row["protein2"]])
                    self.edges[uniprot[row["protein1"]], uniprot[row["protein2"]]][
                        "STRING"
                    ] = (row["combined_score"] / 1000)

                yield (
                    uniprot[row["protein1"]],
                    uniprot[row["protein2"]],
                    self.edges[uniprot[row["protein1"]], uniprot[row["protein2"]]][
                        "STRING"
                    ],
                )

    def add_interactions_from_reactome(
        self, interaction_detection_methods=[], interaction_types=[], reactome_score=0.0
    ):
        for row in download.iterate_tabular_txt(
            data.REACTOME,
            delimiter="\t",
            header=0,
            usecols=[
                "#ID(s) interactor A",
                "ID(s) interactor B",
                "Alt. ID(s) interactor A",
                "Alt. ID(s) interactor B",
                "Taxid interactor A",
                "Taxid interactor B",
                "Interaction detection method(s)",
                "Interaction type(s)",
                "Confidence value(s)",
            ],
        ):

            if not (
                interactor_a := mitab.get_id_from_namespace(
                    row["#ID(s) interactor A"], "uniprotkb"
                )
            ):
                if not (
                    interactor_a := mitab.get_id_from_namespace(
                        row["Alt. ID(s) interactor A"], "uniprotkb"
                    )
                ):
                    continue
            if interactor_a not in self:
                continue

            if not (
                interactor_b := mitab.get_id_from_namespace(
                    row["ID(s) interactor B"], "uniprotkb"
                )
            ):
                if not (
                    interactor_b := mitab.get_id_from_namespace(
                        row["Alt. ID(s) interactor B"], "uniprotkb"
                    )
                ):
                    continue
            if interactor_b not in self:
                continue

            if interactor_a == interactor_b:
                continue

            if not (
                mitab.namespace_has_id(row["Taxid interactor A"], "taxid", "9606")
                and mitab.namespace_has_id(row["Taxid interactor B"], "taxid", "9606")
            ):
                continue

            if not (
                not interaction_detection_methods
                or mitab.namespace_has_any_term_from(
                    row["Interaction detection method(s)"],
                    "psi-mi",
                    interaction_detection_methods,
                )
            ):
                continue

            if not (
                not interaction_types
                or mitab.namespace_has_any_term_from(
                    row["Interaction type(s)"], "psi-mi", interaction_types
                )
            ):
                continue

            if score := mitab.get_id_from_namespace(
                row["Confidence value(s)"], "reactome-score"
            ):
                score = float(score)
                if score < reactome_score:
                    continue
            else:
                continue

            if self.has_edge(interactor_a, interactor_b):
                self.edges[interactor_a, interactor_b]["Reactome"] = max(
                    score, self.edges[interactor_a, interactor_b].get("Reactome", 0.0)
                )
            else:
                self.add_edge(interactor_a, interactor_b)
                self.edges[interactor_a, interactor_b]["Reactome"] = score

            yield (
                interactor_a,
                interactor_b,
                self.edges[interactor_a, interactor_b]["Reactome"],
            )

    def get_modules(self, min_size=3, max_size=100, merge_sizes=max, weight=None):
        communities = [list(self.nodes)]
        while merge_sizes(len(community) for community in communities) > max_size:
            max_indices = [
                communities.index(community)
                for community in communities
                if len(community) == max(len(community) for community in communities)
            ]
            subdivision = False
            for i in range(len(max_indices)):
                subdivided_community = list(
                    modularization.louvain(
                        self.subgraph(communities[max_indices[i]]), weight=weight
                    )
                )
                if len(subdivided_community) > 1:
                    communities[
                        max_indices[i] : max_indices[i] + 1
                    ] = subdivided_community

                    max_indices[i + 1 :] = [
                        index + len(subdivided_community) - 1
                        for index in max_indices[i + 1 :]
                    ]
                    subdivision = True
            if not subdivision:
                break
        return [community for community in communities if len(community) >= min_size]

    def get_proteins(self, time, ptm):
        proteins = []
        for protein in self:
            sites = [
                self.nodes[protein][change]
                for change in self.nodes[protein]
                if len(change.split(" ")) == 3
                and change.split(" ")[0].isnumeric()
                and change.split(" ")[0] == str(time)
                and change.split(" ")[1] == ptm
            ]
            if sites:
                proteins.append(protein)
        return proteins

    def get_proteins_by_change(
        self, time, ptm, change_filter, merge_sites=lambda sites: max(sites, key=abs)
    ):
        proteins = []
        for protein in self:
            sites = [
                self.nodes[protein][change]
                for change in self.nodes[protein]
                if len(change.split(" ")) == 3
                and change.split(" ")[0].isnumeric()
                and change.split(" ")[0] == str(time)
                and change.split(" ")[1] == ptm
            ]
            if sites and change_filter(merge_sites(sites)):
                proteins.append(protein)
        return proteins

    def get_change_enrichment(
        self,
        p=0.05,
        mid_range_thresholds=(-1.0, 1.0),
        get_mid_range=lambda time, ptm, merge_sites, mid_range_thresholds: mid_range_thresholds,
        merge_sites=lambda sites: max(sites, key=abs),
        min_size=3,
        max_size=100,
        weight="weight",
    ):
        modules = self.get_modules(min_size, max_size, weight=weight)
        modules = {i: module for i, module in enumerate(modules)}

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

                mid_range = get_mid_range(
                    time,
                    ptm,
                    merge_sites,
                    mid_range_thresholds,
                )

                n = len(
                    self.get_proteins_by_change(
                        time,
                        ptm,
                        lambda change: change < mid_range[0],
                        merge_sites=merge_sites,
                    )
                ) + len(
                    self.get_proteins_by_change(
                        time,
                        ptm,
                        lambda change: change > mid_range[1],
                        merge_sites=merge_sites,
                    )
                )

                k = [
                    len(
                        self.subgraph(modules[i]).get_proteins_by_change(
                            time,
                            ptm,
                            lambda change: change < mid_range[0],
                            merge_sites=merge_sites,
                        )
                    )
                    + len(
                        self.subgraph(modules[i]).get_proteins_by_change(
                            time,
                            ptm,
                            lambda change: change > mid_range[1],
                            merge_sites=merge_sites,
                        )
                    )
                    for i in modules
                ]

                p_values = {
                    i: scipy.stats.hypergeom.sf(k[i] - 1, M, n, N[i])
                    for i in range(len(modules))
                }

                # Benjamini-Hochberg p-value adjustment (Yekutieli & Benjamini (1999))
                enrichments = sorted(
                    [[p_value, module] for module, p_value in p_values.items()]
                )

                for i in range(len(enrichments)):
                    enrichments[i][0] = min(
                        len(enrichments) * enrichments[i][0] / (i + 1), 1.0
                    )
                    for j in range(i, 0, -1):
                        if enrichments[j - 1][0] <= enrichments[j][0]:
                            break
                        enrichments[j - 1][0] = enrichments[j][0]

                for p_value, module in enrichments:
                    if p_value <= p:
                        if time not in p_adjusted:
                            p_adjusted[time] = {}
                        if ptm not in p_adjusted[time]:
                            p_adjusted[time][ptm] = {}

                        modules_filtered[module] = modules[module]
                        p_adjusted[time][ptm][module] = p_value

        return modules_filtered, p_adjusted
