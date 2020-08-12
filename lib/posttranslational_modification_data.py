import math

import numpy as np
import pandas as pd
import progressbar

from lib.style import CHART_POSITION


def log10_to_log2(x):
    """convert the base of a log from 10 to 2"""
    return math.log10(x) / math.log10(2)


class PosttranslationalModificationData(object):
    """
    class to parse posttranslational modification data
    from one or separate files for each time of measurement
    """

    def __init__(self, ptmtype, configuration):
        # long identifier
        self.ptmtype = ptmtype

        self.configuration = configuration

        # dictionary associating a protein with ratios of deviation
        self.proteins = {}

        # times of measurement
        self.time_steps = self.configuration.get("TIME STEPS")

        # number of sites exhibiting largest deviations to consider
        self.num_sites = self.configuration.get("SITES", 5)

        # short identifier (GraphML node annotation)
        self.id = self.configuration["ID"]

        # color of bar charts in Cytoscape
        self.color = self.configuration.get("COLOR", "#FFFFFF")

        # designation of column reporting UniProt IDs of protein
        self.col_protein = self.configuration.get("COLUMNS", {}).get("PROTEIN")

        # delimiter of sequence positions of modification sites
        self.col_position_sep = self.configuration.get("POSITION DELIMITER",
                                                       ",")

        # delimiter of protein identifiers
        # a peptide may be associated with multiple proteins
        self.col_protein_sep = self.configuration.get("PROTEIN ID DELIMITER")

        # either columns reporting separate times of measurement (single file) or
        # columns reporting separate replicates at one time of measurement (multiple files)
        self.samples = self.configuration.get("SAMPLES")

        # Base of logarithm of ratios of deviation
        self.log_base = self.configuration.get("LOG BASE")

        # designation of column reporting sequence positions
        self.col_position = self.configuration.get("COLUMNS",
                                                   {}).get("POSITION")

        # position of UniProt ID in composite cell value
        self.protein_index = self.configuration.get("COLUMN INDEX PROTEIN ID", 0)

        # Index of header row
        self.header_line = self.configuration.get("HEADER", 0)

        self.significance = float(self.configuration.get("THRESHOLD", "1.0"))

        # minimum number of available replicates
        self.min_replicates = self.configuration.get("REPLICATES", 1)
        if type(self.min_replicates) == str and self.min_replicates.lower(
        ) == "all":
            self.min_replicates = len(self.samples)

        # set function to convert ratio
        if self.log_base is None:
            self.log_conv = math.log2
        elif self.log_base == 2:
            self.log_conv = lambda x: x
        elif self.log_base == 10:
            self.log_conv = log10_to_log2
        else:
            self.log_conv = lambda x: x

        if self.configuration.get("FILE"):
            self.parse_single_file()
        elif self.configuration.get("FILES"):
            self.parse_multiple_files()
        self.merge_site_data()

    def parse_multiple_files(self):
        """
        retrieve posttranslational modification data from multiple files
        each file represents a separate time of measurement
        """
        print("{}:".format(self.ptmtype))
        for file_name, time_step in zip(self.configuration.get("FILES"),
                                        self.configuration.get("TIME STEPS")):
            # set columns to load into a pandas DataFrame
            cols_to_use = [
                self.col_protein,
                self.col_position,
            ]
            # each replicate is reported in a separate column
            cols_to_use.extend(self.samples)
            modification_data = pd.read_excel(io=file_name,
                                              dtype=str,
                                              usecols=cols_to_use,
                                              header=self.header_line)
            bar = progressbar.ProgressBar(max_value=modification_data.shape[0])
            bar.update(0)
            for i, row in modification_data.iterrows():
                # retrieve protein identifier
                protein = row[self.col_protein]

                # proceed with next row if no protein identifier is reported
                if pd.isnull(protein):
                    bar.update(i + 1)
                    continue

                # parse separate replicates
                samples = [
                    float(row[sample])
                    for sample in self.samples
                    if not pd.isnull(row[sample])
                ]

                # if the required number of replicates is not met
                # proceed with next row
                if len(samples) < self.min_replicates:
                    bar.update(i + 1)
                    continue

                # aggregate multiple replicates using their median
                modification = np.median(samples)

                # convert their ratio to a logarithm at a base of 2
                modification = self.log_conv(modification)

                # if no sequence position is found
                # proceed with next row
                if pd.isnull(row[self.col_position]):
                    bar.update(i + 1)
                    continue

                # use first reported position and discard the amino acid
                position = row[self.col_position].split(
                    self.col_position_sep)[0].split("_")[0]

                # add median of measurement to dictionary
                if protein not in self.proteins:
                    self.proteins[protein] = {
                        time_step: {} for time_step in self.time_steps
                    }
                if position not in self.proteins[protein][time_step]:
                    self.proteins[protein][time_step][position] = []

                self.proteins[protein][time_step][position].append(modification)
                bar.update(i + 1)
            bar.finish()

    def parse_single_file(self):
        """
        retrieve posttranslational modification data from a single file
        separate times of measurement are in separate columns
        """
        # set columns to load into a pandas DataFrame
        cols_to_use = [
            self.col_protein,
            self.col_position,
        ]
        cols_to_use.extend(self.samples)
        modification_data = pd.read_excel(io=self.configuration.get("FILE"),
                                          dtype=str,
                                          usecols=cols_to_use,
                                          header=self.header_line)

        print("{}:".format(self.ptmtype))
        bar = progressbar.ProgressBar(max_value=modification_data.shape[0])
        bar.update(0)
        for i, row in modification_data.iterrows():
            # get UniProt ID
            protein = row[self.col_protein].split(
                self.col_protein_sep)[self.protein_index]

            # if no protein identifier is found
            # proceed with next row
            if pd.isnull(protein):
                bar.update(i + 1)
                continue

            # if no sequence position is found
            # proceed with next row
            if pd.isnull(row[self.col_position]):
                bar.update(i + 1)
                continue
            # get sequence position and discard the amino acid
            position = int(row[self.col_position].split(
                self.col_position_sep)[0].split("_")[0])

            # store deviations at individual times of measurement
            modification = {
                time_step: self.log_conv(float(row[sample]))
                for time_step, sample in zip(self.time_steps, self.samples)
                if not pd.isnull(row[sample])
            }
            # associate deviations with UniProt ID
            if protein not in self.proteins:
                self.proteins[protein] = {
                    time_step: {} for time_step in self.time_steps
                }
            for time_step in self.time_steps:
                if position not in self.proteins[protein][time_step]:
                    self.proteins[protein][time_step][position] = []

            for time_step in self.time_steps:
                self.proteins[protein][time_step][position].append(
                    modification[time_step])
            bar.update(i + 1)
        bar.finish()

    def merge_site_data(self):
        # aggregate measurements applying to
        # the same protein, time of measurement and modification site
        # using the median
        for protein in self.proteins:
            for time_step in self.time_steps:
                for modification_site in self.proteins[protein][time_step]:
                    self.proteins[protein][time_step][
                        modification_site] = np.median(
                            self.proteins[protein][time_step]
                            [modification_site])

                # identify sites with largest deviations
                max_pos = sorted(self.proteins[protein][time_step].items(),
                                 key=lambda item: abs(item[1]),
                                 reverse=True)

                # associate a protein with a fixed number of largest deviations at individual sites
                # ordered according to the protein sequence
                # and its single largest deviation
                if max_pos:
                    self.proteins[protein][time_step] = {
                        "SITES": [
                            change
                            for site, change in sorted(max_pos[:self.num_sites])
                        ],
                        "MAX": max([value for site, value in max_pos], key=abs)
                    }
                else:
                    self.proteins[protein][time_step] = {
                        "SITES": [],
                        "MAX": float("nan")
                    }
                # add NaN values should a protein have less modification sites
                # than representable
                self.proteins[protein][time_step]["SITES"].extend([
                    float("nan") for v in
                    range(self.num_sites -
                          len(self.proteins[protein][time_step]["SITES"]))
                ])
