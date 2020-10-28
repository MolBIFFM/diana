import math

import networkx as nx
import numpy as np
import pandas as pd

import ppi.data
import ppi.utilities


class ProteinProteinInteractionNetwork(nx.Graph):
    def __init__(self):
        super(ProteinProteinInteractionNetwork, self).__init__()

    def add_proteins_from_excel(
            self,
            file_name,
            ptm,
            time,
            skiprows=0,
            protein_id_col="Protein",
            protein_id_format=lambda entry: entry,
            position_col="Positions within proteins",
            position_format=lambda entry: entry.split(";")[0].split(".")[0],
            replicates=[
                "Ratio H/L normalized Exp1", "Ratio H/L normalized Exp2",
                "Ratio H/L normalized Exp3"
            ],
            min_num_replicates=1,
            merge_replicates=np.mean,
            convert_measurement=math.log2):
        for _, row in pd.read_excel(file_name,
                                    skiprows=skiprows,
                                    usecols=[protein_id_col, position_col] +
                                    replicates).iterrows():

            if pd.isna(row[protein_id_col]) or pd.isna(row[position_col]):
                continue

            measurements = [
                row[repl] for repl in replicates if not pd.isna(row[repl])
            ]

            if len(measurements) >= min(min_num_replicates, len(replicates)):
                protein_id = str(protein_id_format(str(row[protein_id_col])))
                position = int(position_format(str(row[position_col])))

                if protein_id not in self.nodes:
                    self.add_node(protein_id)

                if position not in self.nodes[protein_id]:
                    self.nodes[protein_id][position] = {}

                if ptm not in self.nodes[protein_id][position]:
                    self.nodes[protein_id][position][ptm] = {}

                self.nodes[protein_id][position][ptm][
                    time] = convert_measurement(merge_replicates(measurements))

    def add_interactions_from_STRING(self,
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
                                     combined_score=0.0):
        protein_ids = {}
        for _, row in ppi.utilities.get_tabular_data(ppi.data.UNIPROT_ID_MAP,
                                                     delimiter="\t",
                                                     usecols=[0, 1, 2]):
            if row[1] == "STRING" and row[0] in self.nodes:
                protein_ids[row[2]] = row[0]

        for _, row in ppi.utilities.get_tabular_data(ppi.data.STRING_ID_MAP,
                                                     usecols=[1, 2]):
            if row[1].split("|")[0] in self.nodes:
                protein_ids[row[2]] = row[1].split("|")[0]

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
                "combined_score": combined_score
            }.items() if threshold
        }

        for _, row in ppi.utilities.get_tabular_data(
                ppi.data.STRING,
                delimiter=" ",
                header=0,
                usecols=["protein1", "protein2"] + list(thresholds.keys())):
            if (protein_ids.get(row["protein1"])
                    and protein_ids.get(row["protein2"])
                    and all(row[column] / 1000 >= thresholds[column]
                            for column in thresholds)):
                self.add_edge(protein_ids[row["protein1"]],
                              protein_ids[row["protein2"]])
