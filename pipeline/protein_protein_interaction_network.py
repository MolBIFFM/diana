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
            header=0,
            protein_id_col="Protein",
            protein_id_format=lambda entry: entry,
            position_col="Positions within proteins",
            position_format=lambda entry: entry.split(";")[0].split(".")[0],
            replicates=[
                "Ratio H/L normalized Exp1", "Ratio H/L normalized Exp2",
                "Ratio H/L normalized Exp3"
            ],
            min_num_replicates=2,
            merge_replicates=statistics.mean,
            convert_measurement=math.log2,
            num_sites=5):
        proteins = {}
        for _, row in pd.read_excel(
                file_name,
                header=header,
                usecols=[protein_id_col, position_col] + replicates,
                dtype={
                    protein_id_col: str,
                    position_col: str,
                    **{replicate: float
                       for replicate in replicates}
                }).iterrows():

            if pd.isna(row[protein_id_col]) or pd.isna(row[position_col]):
                continue

            measurements = [
                row[repl] for repl in replicates if not pd.isna(row[repl])
            ]

            if len(measurements) >= min(min_num_replicates, len(replicates)):
                protein_id = str(protein_id_format(row[protein_id_col]))
                position = int(position_format(row[position_col]))

                if protein_id not in proteins:
                    proteins[protein_id] = []
                bisect.insort(
                    proteins[protein_id],
                    (position,
                     convert_measurement(merge_replicates(measurements))))

        for protein in proteins:
            self.add_node(protein)
            if len(proteins[protein]) > num_sites:
                proteins[protein] = sorted(
                    sorted(proteins[protein],
                           key=lambda tp: tp[1],
                           reverse=True)[:num_sites])
            for i in range(len(proteins[protein])):
                self.nodes[protein][" ".join(
                    ["site", ptm, str(time),
                     str(i + 1)])] = proteins[protein][i][0]
                self.nodes[protein][" ".join(
                    ["change", ptm, str(time),
                     str(i + 1)])] = proteins[protein][i][1]

    def post_translational_modifications(self):
        return tuple(
            sorted(
                set(
                    attr.split(" ")[1] for protein in self
                    for attr in self.nodes[protein]
                    if len(attr.split(" ")) == 4
                    and attr.split(" ")[0] == "change")))

    def times(self):
        return tuple(
            sorted(
                set(
                    int(attr.split(" ")[2]) for protein in self
                    for attr in self.nodes[protein]
                    if len(attr.split(" ")) == 4
                    and attr.split(" ")[0] == "change")))

    def cytoscape_style_shape(self):
        modifications = {
            ("P", "U"): "both",
            ("P", ): "phosphorylation",
            ("U", ): "ubiquitination",
            (): "none"
        }

        for time in self.times():
            for protein in self:
                self.nodes[protein][
                    "post-translational modification {}".format(
                        time)] = modifications[tuple(
                            sorted(
                                set(
                                    attr.split(" ")[1]
                                    for attr in self.nodes[protein]
                                    if len(attr.split(" ")) == 4
                                    and attr.split(" ")[0] == "change")))]

    def cytoscape_style_color(self,
                              merge_changes=statistics.mean,
                              threshold=1.0):
        post_translational_modifications = self.post_translational_modifications(
        )
        for time in self.times():
            for protein in self:
                ptm = {}
                for post_translational_modification in post_translational_modifications:
                    changes = [
                        self.nodes[protein][attr]
                        for attr in self.nodes[protein]
                        if len(attr.split(" ")) == 4
                        and attr.split(" ")[0] == "change" and attr.split(
                            " ")[1] == post_translational_modification
                        and attr.split(" ")[2] == str(time)
                    ]
                    if changes:
                        ptm[post_translational_modification] = merge_changes(
                            changes)

                if all(change > 0.0 for change in ptm.values()):
                    if any(change >= threshold for change in ptm.values()):
                        self.nodes[protein]["change {}".format(time)] = "up"
                    else:
                        self.nodes[protein]["change {}".format(
                            time)] = "mid up"
                elif all(change < 0.0 for change in ptm.values()):
                    if any(change <= -threshold for change in ptm.values()):
                        self.nodes[protein]["change {}".format(time)] = "down"
                    else:
                        self.nodes[protein]["change {}".format(
                            time)] = "mid down"
                else:
                    if ptm:
                        self.nodes[protein]["change {}".format(
                            time)] = " ".join([
                                "{} up".format(post_translational_modification)
                                if ptm[post_translational_modification] > 0.0
                                else "{} down".format(
                                    post_translational_modification)
                                for post_translational_modification in ptm
                            ])
                    else:
                        self.nodes[protein]["change {}".format(time)] = "mid"

    def add_interactions_from_BioGRID(
        self,
        experimental_system=[
            "Affinity Capture-Luminescence", "Affinity Capture-MS",
            "Affinity Capture-RNA", "Affinity Capture-Western",
            "Biochemical Activity", "Co-crystal Structure", "Co-purification",
            "FRET", "PCA", "Two-hybrid"
        ]):
        uniprot = {}
        for row in fetch.read_tabular_data(
                protein_protein_interaction_data.UNIPROT_ID_MAP,
                delimiter="\t",
                usecols=[0, 1, 2]):
            if row[1] == "BioGRID" and row[0] in self.nodes:
                uniprot[int(row[2])] = row[0]

        for row in fetch.read_tabular_data(
                protein_protein_interaction_data.BIOGRID_ID_MAP_ARCHIVE,
                zip_file=protein_protein_interaction_data.BIOGRID_ID_MAP_FILE,
                delimiter="\t",
                header=20,
                usecols=[
                    "BIOGRID_ID", "IDENTIFIER_VALUE", "IDENTIFIER_TYPE",
                    "ORGANISM_OFFICIAL_NAME"
                ]):
            if (row["IDENTIFIER_TYPE"]
                    in ("UNIPROT-ACCESSION", "UNIPROT-ISOFORM")
                    and row["ORGANISM_OFFICIAL_NAME"] == "Homo sapiens"
                    and row["IDENTIFIER_VALUE"] in self.nodes):
                uniprot[int(row["BIOGRID_ID"])] = row["IDENTIFIER_VALUE"]

        for row in fetch.read_tabular_data(
                protein_protein_interaction_data.BIOGRID_TAB3_ARCHIVE,
                zip_file=protein_protein_interaction_data.BIOGRID_TAB3_FILE,
                delimiter="\t",
                header=0,
                usecols=[
                    "BioGRID ID Interactor A", "BioGRID ID Interactor B",
                    "Experimental System", "Experimental System Type",
                    "Organism ID Interactor A", "Organism ID Interactor B",
                    "Throughput"
                ]):
            if (uniprot.get(row["BioGRID ID Interactor A"])
                    and uniprot.get(row["BioGRID ID Interactor B"])
                    and row["BioGRID ID Interactor A"] !=
                    row["BioGRID ID Interactor B"]
                    and row["Experimental System"] in experimental_system
                    and row["Experimental System Type"] == "physical"
                    and row["Organism ID Interactor A"] == 9606
                    and row["Organism ID Interactor B"] == 9606):
                self.add_edge(uniprot[row["BioGRID ID Interactor A"]],
                              uniprot[row["BioGRID ID Interactor B"]])
                self.edges[
                    uniprot[row["BioGRID ID Interactor A"]],
                    uniprot[row["BioGRID ID Interactor B"]]]["BioGRID"] = 1.0

    def add_interactions_from_BioGRID_MITAB(
        self,
        interaction_detection_methods=[
            "affinity chromatography technology", "two hybrid", "biochemical",
            "pull down", "enzymatic study", "bio id", "x-ray crystallography",
            "fluorescent resonance energy transfer",
            "protein complementation assay"
        ],
        interaction_types=[
            "physical association", "direct interaction", "association"
        ]):
        uniprot = {}
        for row in fetch.read_tabular_data(
                protein_protein_interaction_data.UNIPROT_ID_MAP,
                delimiter="\t",
                usecols=[0, 1, 2]):
            if row[1] == "BioGRID" and row[0] in self.nodes:
                uniprot[row[2]] = row[0]

        for row in fetch.read_tabular_data(
                protein_protein_interaction_data.BIOGRID_ID_MAP_ARCHIVE,
                zip_file=protein_protein_interaction_data.BIOGRID_ID_MAP_FILE,
                delimiter="\t",
                header=20,
                usecols=[
                    "BIOGRID_ID", "IDENTIFIER_VALUE", "IDENTIFIER_TYPE",
                    "ORGANISM_OFFICIAL_NAME"
                ]):
            if (row["IDENTIFIER_TYPE"]
                    in ("UNIPROT-ACCESSION", "UNIPROT-ISOFORM")
                    and row["ORGANISM_OFFICIAL_NAME"] == "Homo sapiens"
                    and row["IDENTIFIER_VALUE"] in self.nodes):
                uniprot[int(row["BIOGRID_ID"])] = row["IDENTIFIER_VALUE"]

        for row in fetch.read_tabular_data(
                protein_protein_interaction_data.BIOGRID_MITAB_ARCHIVE,
                zip_file=protein_protein_interaction_data.BIOGRID_MITAB_FILE,
                delimiter="\t",
                header=0,
                usecols=[
                    "#ID Interactor A", "ID Interactor B",
                    "Alt IDs Interactor A", "Alt IDs Interactor B",
                    "Taxid Interactor A", "Taxid Interactor B",
                    "Interaction Detection Method", "Interaction Types"
                ]):
            if not (interactor_a := mitab.get_id_from_ns(
                    row["#ID Interactor A"], "biogrid")):
                if not (interactor_a := mitab.get_id_from_ns(
                        row["Alt IDs Interactor A"], "biogrid")):
                    continue
            if not uniprot.get(interactor_a):
                continue

            if not (interactor_b := mitab.get_id_from_ns(
                    row["ID Interactor B"], "biogrid")):
                if not (interactor_b := mitab.get_id_from_ns(
                        row["Alt IDs Interactor B"], "biogrid")):
                    continue
            if not uniprot.get(interactor_b):
                continue

            if interactor_a == interactor_b:
                continue

            if not (mitab.ns_has_id(row["Taxid Interactor A"], "taxid", "9606")
                    and mitab.ns_has_id(row["Taxid Interactor B"], "taxid",
                                        "9606")):
                continue

            if not mitab.ns_has_any_term(row["Interaction Detection Method"],
                                         "psi-mi",
                                         interaction_detection_methods):
                continue

            if not mitab.ns_has_any_term(row["Interaction Types"], "psi-mi",
                                         interaction_types):
                continue

            self.add_edge(uniprot[interactor_a], uniprot[interactor_b])
            self.edges[interactor_a, interactor_b]["BioGRID"] = 1.0

    def add_interactions_from_IntAct(
            self,
            interaction_detection_methods=[
                "affinity chromatography technology", "two hybrid",
                "biochemical", "pull down", "enzymatic study", "bio id",
                "x-ray crystallography",
                "fluorescent resonance energy transfer",
                "protein complementation assay"
            ],
            interaction_types=[
                "physical association", "direct interaction", "association"
            ],
            miscore=0.27):
        for row in fetch.read_tabular_data(
                protein_protein_interaction_data.INTACT_ARCHIVE,
                zip_file=protein_protein_interaction_data.INTACT_FILE,
                delimiter="\t",
                header=0,
                usecols=[
                    "#ID(s) interactor A", "ID(s) interactor B",
                    "Alt. ID(s) interactor A", "Alt. ID(s) interactor B",
                    "Taxid interactor A", "Taxid interactor B",
                    "Interaction detection method(s)", "Interaction type(s)",
                    "Confidence value(s)"
                ]):

            if not (interactor_a := mitab.get_id_from_ns(
                    row["#ID(s) interactor A"], "uniprotkb")):
                if not (interactor_a := mitab.get_id_from_ns(
                        row["Alt. ID(s) interactor A"], "uniprotkb")):
                    continue
            if interactor_a not in self:
                continue

            if not (interactor_b := mitab.get_id_from_ns(
                    row["ID(s) interactor B"], "uniprotkb")):
                if not (interactor_b := mitab.get_id_from_ns(
                        row["Alt. ID(s) interactor B"], "uniprotkb")):
                    continue
            if interactor_b not in self:
                continue

            if interactor_a == interactor_b:
                continue

            if not (mitab.ns_has_id(row["Taxid interactor A"], "taxid", "9606")
                    and mitab.ns_has_id(row["Taxid interactor B"], "taxid",
                                        "9606")):
                continue

            if not mitab.ns_has_any_term(
                    row["Interaction detection method(s)"], "psi-mi",
                    interaction_detection_methods):
                continue

            if not mitab.ns_has_any_term(row["Interaction type(s)"], "psi-mi",
                                         interaction_types):
                continue

            if score := mitab.get_id_from_ns(row["Confidence value(s)"],
                                             "intact-miscore"):
                score = float(score)
                if score < miscore:
                    continue
            else:
                continue

            self.add_edge(interactor_a, interactor_b)
            if self.has_edge(interactor_a, interactor_b):
                self.edges[interactor_a, interactor_b]["IntAct"] = max(
                    score, self.edges[interactor_a,
                                      interactor_b].get("IntAct", 0.0))
            else:
                self.edges[interactor_a, interactor_b]["IntAct"] = score

    def add_interactions_from_STRING(self,
                                     neighborhood=0.0,
                                     neighborhood_transferred=0.0,
                                     fusion=0.0,
                                     cooccurence=0.0,
                                     homology=0.0,
                                     coexpression=0.0,
                                     coexpression_transferred=0.0,
                                     experiments=0.7,
                                     experiments_transferred=0.0,
                                     database=0.0,
                                     database_transferred=0.0,
                                     textmining=0.0,
                                     textmining_transferred=0.0,
                                     combined_score=0.7):

        uniprot = {}
        for row in fetch.read_tabular_data(
                protein_protein_interaction_data.UNIPROT_ID_MAP,
                delimiter="\t",
                usecols=[0, 1, 2]):
            if row[1] == "STRING" and row[0] in self.nodes:
                uniprot[row[2]] = row[0]

        for row in fetch.read_tabular_data(
                protein_protein_interaction_data.STRING_ID_MAP, usecols=[1,
                                                                         2]):
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
                "textmining_transferred": textmining_transferred
            }.items() if threshold
        }
        thresholds["combined_score"] = combined_score

        for row in fetch.read_tabular_data(
                protein_protein_interaction_data.STRING,
                delimiter=" ",
                header=0,
                usecols=["protein1", "protein2"] + list(thresholds.keys())):
            if (uniprot.get(row["protein1"]) and uniprot.get(row["protein2"])
                    and row["protein1"] != row["protein2"]
                    and all(row[column] / 1000 >= thresholds[column]
                            for column in thresholds)):
                self.add_edge(uniprot[row["protein1"]],
                              uniprot[row["protein2"]])
                if self.has_edge(uniprot[row["protein1"]],
                                 uniprot[row["protein2"]]):
                    self.edges[uniprot[row["protein1"]],
                               uniprot[row["protein2"]]]["STRING"] = max(
                                   row["combined_score"] / 1000,
                                   self.edges[uniprot[row["protein1"]],
                                              uniprot[row["protein2"]]].get(
                                                  "STRING", 0.0))
                else:
                    self.edges[uniprot[row["protein1"]],
                               uniprot[row["protein2"]]][
                                   "STRING"] = row["combined_score"] / 1000

    def remove_isolates(self):
        self.remove_nodes_from(list(nx.isolates(self)))

    def export_as_graphml(self, file_name):
        nx.write_graphml_xml(self, file_name)

    def export_as_cyjs(self, file_name):
        with open(file_name, "w") as file:
            json.dump(nx.readwrite.json_graph.cytoscape_data(self),
                      file,
                      indent=2)
