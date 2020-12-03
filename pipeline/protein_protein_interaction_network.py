import math
import statistics
import bisect
import json

import networkx as nx
import pandas as pd

from pipeline.configuration import protein_protein_interaction_data
from pipeline.utilities import fetch
from pipeline.utilities import mitab


class ProteinProteinInteractionNetwork(nx.Graph):
    def __init__(self):
        super(ProteinProteinInteractionNetwork, self).__init__()

    def add_proteins_from_excel(
        self,
        file_name,
        ptm,
        time,
        protein_id_col,
        position_col,
        replicates,
        protein_id_format,
        position_format,
        header=0,
        num_sites=1,
        num_replicates=1,
        merge_replicates=statistics.mean,
        convert_measurement=math.log2,
    ):
        proteins = {}
        for _, row in pd.read_excel(
            file_name,
            header=header,
            usecols=[protein_id_col, position_col] + replicates,
            dtype={
                protein_id_col: str,
                position_col: str,
                **{replicate: float for replicate in replicates},
            },
        ).iterrows():

            if pd.isna(row[protein_id_col]) or pd.isna(row[position_col]):
                continue

            measurements = [row[repl] for repl in replicates if not pd.isna(row[repl])]

            if len(measurements) >= min(num_replicates, len(replicates)):
                protein_id = str(protein_id_format(row[protein_id_col]))
                position = int(position_format(row[position_col]))

                if protein_id not in proteins:
                    proteins[protein_id] = []
                bisect.insort(
                    proteins[protein_id],
                    (position, convert_measurement(merge_replicates(measurements))),
                )

        for protein in proteins:
            self.add_node(protein)
            if len(proteins[protein]) > num_sites:
                proteins[protein] = sorted(
                    sorted(proteins[protein], key=lambda tp: tp[1], reverse=True)[
                        :num_sites
                    ]
                )
            for i in range(len(proteins[protein])):
                self.nodes[protein][
                    "{} {} {}".format(str(time), ptm, str(i + 1))
                ] = proteins[protein][i][1]

            yield (
                protein,
                tuple([proteins[protein][i][1] for i in range(len(proteins[protein]))]),
            )

    def get_times(self):
        return tuple(
            sorted(
                set(
                    int(attr.split(" ")[0])
                    for protein in self
                    for attr in self.nodes[protein]
                    if len(attr.split(" ")) == 3
                )
            )
        )

    def get_post_translational_modifications(self):
        return {
            time: tuple(
                sorted(
                    set(
                        attr.split(" ")[1]
                        for protein in self
                        for attr in self.nodes[protein]
                        if len(attr.split(" ")) == 3 and attr.split(" ")[0] == str(time)
                    )
                )
            )
            for time in self.get_times()
        }

    def get_sites(self):
        return {
            time: {
                ptm: max(
                    int(attr.split(" ")[2])
                    for protein in self
                    for attr in self.nodes[protein]
                    if len(attr.split(" ")) == 3
                    and attr.split(" ")[0] == str(time)
                    and attr.split(" ")[1] == ptm
                )
                for ptm in self.get_post_translational_modifications()[time]
            }
            for time in self.get_times()
        }

    def set_ptm_data_column(self):
        for time in self.get_times():
            for protein in self:
                self.nodes[protein]["PTM {}".format(time)] = " ".join(
                    tuple(
                        sorted(
                            set(
                                attr.split(" ")[1]
                                for attr in self.nodes[protein]
                                if len(attr.split(" ")) == 3
                                and attr.split(" ")[0] == str(time)
                            )
                        )
                    )
                )

    def set_trend_data_column(
        self, merge_trends=statistics.mean, mid_range=(-1.0, 1.0)
    ):
        modifications = self.get_post_translational_modifications()
        for time in self.get_times():
            for protein in self:
                ptm = {}
                for post_translational_modification in modifications[time]:
                    trends = [
                        self.nodes[protein][attr]
                        for attr in self.nodes[protein]
                        if len(attr.split(" ")) == 3
                        and attr.split(" ")[0] == str(time)
                        and attr.split(" ")[1] == post_translational_modification
                    ]
                    if trends:
                        ptm[post_translational_modification] = merge_trends(trends)

                if ptm:
                    if all(trend > 0.0 for trend in ptm.values()):
                        if any(trend >= mid_range[1] for trend in ptm.values()):
                            self.nodes[protein]["trend {}".format(time)] = "up"
                        else:
                            self.nodes[protein]["trend {}".format(time)] = "mid up"
                    elif all(trend < 0.0 for trend in ptm.values()):
                        if any(trend <= mid_range[0] for trend in ptm.values()):
                            self.nodes[protein]["trend {}".format(time)] = "down"
                        else:
                            self.nodes[protein]["trend {}".format(time)] = "mid down"
                    else:
                        self.nodes[protein]["trend {}".format(time)] = " ".join(
                            sorted(
                                [
                                    "{} up".format(post_translational_modification)
                                    if ptm[post_translational_modification] > 0.0
                                    else "{} down".format(
                                        post_translational_modification
                                    )
                                    for post_translational_modification in ptm
                                ]
                            )
                        )
                else:
                    self.nodes[protein]["trend {}".format(time)] = ""

    def add_interactions_from_BioGRID(self, experimental_system=[]):
        uniprot = {}
        for row in fetch.read_tabular_data(
            protein_protein_interaction_data.UNIPROT_ID_MAP,
            delimiter="\t",
            usecols=[0, 1, 2],
        ):
            if row[1] == "BioGRID" and row[0] in self.nodes:
                uniprot[int(row[2])] = row[0]

        for row in fetch.read_tabular_data(
            protein_protein_interaction_data.BIOGRID_ID_MAP_ARCHIVE,
            zip_file=protein_protein_interaction_data.BIOGRID_ID_MAP_FILE,
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

        for row in fetch.read_tabular_data(
            protein_protein_interaction_data.BIOGRID_ARCHIVE,
            zip_file=protein_protein_interaction_data.BIOGRID_FILE,
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
        for row in fetch.read_tabular_data(
            protein_protein_interaction_data.INTACT_ARCHIVE,
            zip_file=protein_protein_interaction_data.INTACT_FILE,
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
    ):

        uniprot = {}
        for row in fetch.read_tabular_data(
            protein_protein_interaction_data.UNIPROT_ID_MAP,
            delimiter="\t",
            usecols=[0, 1, 2],
        ):
            if row[1] == "STRING" and row[0] in self.nodes:
                uniprot[row[2]] = row[0]

        for row in fetch.read_tabular_data(
            protein_protein_interaction_data.STRING_ID_MAP, usecols=[1, 2]
        ):
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

        for row in fetch.read_tabular_data(
            protein_protein_interaction_data.STRING,
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

    def remove_isolates(self):
        self.remove_nodes_from(list(nx.isolates(self)))

    def export_as_graphml(self, file_name):
        nx.write_graphml_xml(self, file_name)

    def export_as_cyjs(self, file_name):
        with open(file_name, "w") as file:
            json.dump(nx.readwrite.json_graph.cytoscape_data(self), file, indent=2)
