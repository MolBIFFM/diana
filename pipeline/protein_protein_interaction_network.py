import bisect
import itertools
import json
import math
import statistics

import networkx as nx
import scipy.stats
import pandas as pd

from pipeline.configuration import data
from pipeline.utilities import fetch, mitab


class ProteinProteinInteractionNetwork(nx.Graph):
    def __init__(self):
        super(ProteinProteinInteractionNetwork, self).__init__()

    def add_proteins_from_excel(
        self,
        file_name,
        ptm,
        time,
        protein_accession_col,
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
            usecols=[protein_accession_col, position_col] + replicates,
            dtype={
                protein_accession_col: str,
                position_col: str,
                **{replicate: float for replicate in replicates},
            },
        ).iterrows():

            if pd.isna(row[protein_accession_col]) or pd.isna(row[position_col]):
                continue

            measurements = [row[repl] for repl in replicates if not pd.isna(row[repl])]

            if len(measurements) >= min(num_replicates, len(replicates)):
                protein_accessions = [
                    str(protein_accession)
                    for protein_accession in protein_accession_format(
                        row[protein_accession_col]
                    )
                ]
                positions = [
                    int(position) for position in position_format(row[position_col])
                ]

                for protein_accession, position in itertools.zip_longest(
                    protein_accessions, positions, fillvalue=0
                ):
                    if len(protein_accession.split("-")) > 1:
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

        reviewed_proteins, primary_accession, gene_name = {}, {}, {}
        accession_in_proteins = None
        for line in fetch.iterate_data(data.UNIPROT_SWISS_PROT):
            if line.split("   ")[0] == "AC":
                accessions = line.split("   ")[1].rstrip(";").split("; ")
                for i, accession in enumerate(accessions):
                    if accession in proteins:
                        reviewed_proteins[accession] = proteins[accession]
                        if i > 0:
                            primary_accession[accession] = accessions[0]
                        accession_in_proteins = accession
            elif line.split("   ")[0] == "GN" and accession_in_proteins:
                for entry in line.split("   ")[1].rstrip(";").split("; "):
                    if entry.split(" ")[0].split("=")[0] == "Name":
                        gene_name[accession_in_proteins] = entry.split(" ")[0].split(
                            "="
                        )[1]
                        if accession_in_proteins in primary_accession:
                            gene_name[
                                primary_accession[accession_in_proteins]
                            ] = gene_name[accession_in_proteins]
            elif line == "//":
                accession_in_proteins = None

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
                if protein in gene_name:
                    gene_name[primary_accession[protein]] = gene_name.pop(protein)

        for protein in proteins:
            for isoform in proteins[protein]:
                if isoform != "1":
                    self.add_node("-".join([protein, isoform]))
                    self.nodes["-".join([protein, isoform])][
                        "gene name"
                    ] = gene_name.get(protein, "")
                else:
                    self.add_node(protein)
                    self.nodes[protein]["gene name"] = gene_name.get(protein, "")

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
                    gene_name.get(protein, ""),
                    "-".join([protein, isoform]) if isoform != "1" else protein,
                    [
                        proteins[protein][isoform][i][1]
                        for i in range(len(proteins[protein][isoform]))
                    ],
                ]

    def get_times(self):
        return sorted(
            set(
                int(attribute.split(" ")[0])
                for protein in self
                for attribute in self.nodes[protein]
                if len(attribute.split(" ")) == 3
            )
        )

    def get_post_translational_modifications(self, time):
        return sorted(
            set(
                attribute.split(" ")[1]
                for protein in self
                for attribute in self.nodes[protein]
                if len(attribute.split(" ")) == 3
                and attribute.split(" ")[0] == str(time)
            )
        )

    def get_sites(self, time, ptm):
        return max(
            int(attribute.split(" ")[2])
            for protein in self
            for attribute in self.nodes[protein]
            if len(attribute.split(" ")) == 3
            and attribute.split(" ")[0] == str(time)
            and attribute.split(" ")[1] == ptm
        )

    def set_post_translational_modification_data_column(self):
        for time in self.get_times():
            for protein in self:
                self.nodes[protein]["PTM {}".format(time)] = " ".join(
                    sorted(
                        set(
                            attribute.split(" ")[1]
                            for attribute in self.nodes[protein]
                            if len(attribute.split(" ")) == 3
                            and attribute.split(" ")[0] == str(time)
                        )
                    )
                )

    def get_changes(self, time, ptm, merge_sites=statistics.mean):
        changes = []
        for protein in self:
            sites = [
                self.nodes[protein][attribute]
                for attribute in self.nodes[protein]
                if len(attribute.split(" ")) == 3
                and attribute.split(" ")[0] == str(time)
                and attribute.split(" ")[1] == ptm
            ]

            if sites:
                changes.append(merge_sites(sites))
            else:
                changes.append(0.0)

        return changes

    def get_range_p(self, time, ptm, merge_sites=statistics.mean, p=0.05):
        changes = self.get_changes(time, ptm, merge_sites)

        change_distribution = scipy.stats.norm(
            loc=0.0, scale=scipy.stats.norm.fit(changes, floc=0.0)[1]
        )

        return (
            change_distribution.ppf(0.5 * p),
            change_distribution.ppf(1.0 - 0.5 * p),
        )

    def set_change_data_column_p(self, merge_sites=statistics.mean, p=0.05):
        times = self.get_times()
        post_translational_modifications = {
            time: self.get_post_translational_modifications(time) for time in times
        }

        mid_range = {
            time: {
                post_translational_modification: self.get_range_p(
                    time, post_translational_modification, merge_sites, p
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
                        self.nodes[protein][attribute]
                        for attribute in self.nodes[protein]
                        if len(attribute.split(" ")) == 3
                        and attribute.split(" ")[0] == str(time)
                        and attribute.split(" ")[1] == post_translational_modification
                    ]
                    if sites:
                        ptms[post_translational_modification] = merge_sites(sites)

                if ptms:
                    if all(change > 0.0 for change in ptms.values()):
                        if any(
                            change
                            >= mid_range[time][post_translational_modification][1]
                            for post_translational_modification, change in ptms.items()
                        ):
                            self.nodes[protein]["change {}".format(time)] = "up"
                        else:
                            self.nodes[protein]["change {}".format(time)] = "mid up"
                    elif all(change < 0.0 for change in ptms.values()):
                        if any(
                            change <= mid_range[time][ptm][0]
                            for ptm, change in ptms.items()
                        ):
                            self.nodes[protein]["change {}".format(time)] = "down"
                        else:
                            self.nodes[protein]["change {}".format(time)] = "mid down"
                    else:
                        self.nodes[protein]["change {}".format(time)] = " ".join(
                            sorted(
                                [
                                    ""
                                    if change == 0.0
                                    else "{} up".format(ptm)
                                    if change > 0.0
                                    else "{} down".format(ptm)
                                    for ptm, change in ptms.items()
                                ]
                            )
                        )
                else:
                    self.nodes[protein]["change {}".format(time)] = ""

    def set_change_data_column_abs(
        self, merge_sites=statistics.mean, absolute_change=1.0
    ):
        times = self.get_times()
        post_translational_modifications = {
            time: self.get_post_translational_modifications(time) for time in times
        }

        for time in times:
            for protein in self:
                ptms = {}
                for post_translational_modification in post_translational_modifications[
                    time
                ]:
                    sites = [
                        self.nodes[protein][attribute]
                        for attribute in self.nodes[protein]
                        if len(attribute.split(" ")) == 3
                        and attribute.split(" ")[0] == str(time)
                        and attribute.split(" ")[1] == post_translational_modification
                    ]
                    if sites:
                        ptms[post_translational_modification] = merge_sites(sites)

                if ptms:
                    if all(change > 0.0 for change in ptms.values()):
                        if any(
                            change >= absolute_change for ptm, change in ptms.items()
                        ):
                            self.nodes[protein]["change {}".format(time)] = "up"
                        else:
                            self.nodes[protein]["change {}".format(time)] = "mid up"
                    elif all(change < 0.0 for change in ptms.values()):
                        if any(
                            change <= -absolute_change for ptm, change in ptms.items()
                        ):
                            self.nodes[protein]["change {}".format(time)] = "down"
                        else:
                            self.nodes[protein]["change {}".format(time)] = "mid down"
                    else:
                        self.nodes[protein]["change {}".format(time)] = " ".join(
                            sorted(
                                [
                                    ""
                                    if change == 0.0
                                    else "{} up".format(ptm)
                                    if change > 0.0
                                    else "{} down".format(ptm)
                                    for ptm, change in ptms.items()
                                ]
                            )
                        )
                else:
                    self.nodes[protein]["change {}".format(time)] = ""

    def set_edge_weights(self, weight=lambda confidence_scores: 1.0):
        for edge in self.edges:
            self.edges[edge]["weight"] = weight(
                {
                    database: self.edges[edge][database]
                    for database in ("BioGRID", "IntAct", "STRING")
                    if database in self.edges[edge]
                }
            )

    def add_interactions_from_BioGRID(
        self, experimental_system=[], multi_validated_physical=False
    ):
        uniprot = {}
        for row in fetch.iterate_tabular_data(
            data.UNIPROT_ID_MAP,
            delimiter="\t",
            usecols=[0, 1, 2],
        ):
            if row[1] == "BioGRID" and row[0] in self.nodes:
                uniprot[int(row[2])] = row[0]

        for row in fetch.iterate_tabular_data(
            data.BIOGRID_ID_MAP_ARCHIVE,
            zip_file=data.BIOGRID_ID_MAP_FILE,
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

        for row in fetch.iterate_tabular_data(
            data.BIOGRID_MV_PHYSICAL_ARCHIVE
            if multi_validated_physical
            else data.BIOGRID_ARCHIVE,
            zip_file=data.BIOGRID_MV_PHYSICAL_FILE
            if multi_validated_physical
            else data.BIOGRID_FILE,
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
                ]["BioGRID"] = 1.0

                yield (
                    uniprot[row["BioGRID ID Interactor A"]],
                    uniprot[row["BioGRID ID Interactor B"]],
                )

    def add_interactions_from_IntAct(
        self, interaction_detection_methods=[], interaction_types=[], mi_score=0.0
    ):
        for row in fetch.iterate_tabular_data(
            data.INTACT_ARCHIVE,
            zip_file=data.INTACT_FILE,
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
                interactor_a := mitab.get_id_from_ns(
                    row["#ID(s) interactor A"], "uniprotkb"
                )
            ):
                if not (
                    interactor_a := mitab.get_id_from_ns(
                        row["Alt. ID(s) interactor A"], "uniprotkb"
                    )
                ):
                    continue
            if interactor_a not in self:
                continue

            if not (
                interactor_b := mitab.get_id_from_ns(
                    row["ID(s) interactor B"], "uniprotkb"
                )
            ):
                if not (
                    interactor_b := mitab.get_id_from_ns(
                        row["Alt. ID(s) interactor B"], "uniprotkb"
                    )
                ):
                    continue
            if interactor_b not in self:
                continue

            if interactor_a == interactor_b:
                continue

            if not (
                mitab.ns_has_id(row["Taxid interactor A"], "taxid", "9606")
                and mitab.ns_has_id(row["Taxid interactor B"], "taxid", "9606")
            ):
                continue

            if not (
                not interaction_detection_methods
                or mitab.ns_has_any_term(
                    row["Interaction detection method(s)"],
                    "psi-mi",
                    interaction_detection_methods,
                )
            ):
                continue

            if not (
                not interaction_types
                or mitab.ns_has_any_term(
                    row["Interaction type(s)"], "psi-mi", interaction_types
                )
            ):
                continue

            if score := mitab.get_id_from_ns(
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

    def add_interactions_from_STRING(
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
        for row in fetch.iterate_tabular_data(
            data.UNIPROT_ID_MAP,
            delimiter="\t",
            usecols=[0, 1, 2],
        ):
            if row[1] == "STRING" and row[0] in self.nodes:
                uniprot[row[2]] = row[0]

        for row in fetch.iterate_tabular_data(data.STRING_ID_MAP, usecols=[1, 2]):
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

        for row in fetch.iterate_tabular_data(
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
