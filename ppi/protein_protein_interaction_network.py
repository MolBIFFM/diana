import math

import networkx as nx
import numpy as np
import pandas as pd

import ppi.configuration


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

                if ptm not in self.nodes[protein_id]:
                    self.nodes[protein_id][ptm] = {}

                if position not in self.nodes[protein_id][ptm]:
                    self.nodes[protein_id][ptm][position] = {}

                self.nodes[protein_id][ptm][position][
                    time] = convert_measurement(merge_replicates(measurements))
