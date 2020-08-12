import csv
import ftplib
import gzip
import logging
import math
import os
import pathlib
import re
import zipfile
from urllib.parse import urlparse

import matplotlib
import matplotlib.pylab as pl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import progressbar
import requests
from matplotlib import ticker
from scipy import stats

DPI = 300


def get_psimi(string):
    """Extract PSI-MI ID ("MI:[0-9]{4}") from string"""
    m = re.search(r"(MI:[0-9]{4})", string)
    return m.group(1) if m is not None else None


class InteractionDatabase(object):
    """
    Class to obtain protein-protein interactions 
    from tabular representation of a protein-protein interaction database 
    in MITAB format
    """

    def __init__(self, name, configuration):
        self.logger = logging.getLogger("main")

        # name of the data base
        self.name = name

        # dictionary storing database-specific user configuration
        self.configuration = configuration

        # dictionary associating undirected protein-protein interactions with their confidence scores
        self.interactions = {}

        # specialized container for memory efficient storage of confidence score distribution
        self.scores = self.Scores()

        # set of proteins unsuccessfully mapped to UniProt IDs for logging
        self.unmapped_proteins = set()

        # number of interactions discarded based on score thresholds
        self.discarded_score = 0

        # number of interactions discarded based on UniProt ID mapping
        self.discarded_unmapped = 0

        # number of interactions discarded based on required detection methods or interaction types
        self.discarded_methods_types = 0

        # total number of applicable (directed or redundant) interactions obtained from database
        self.num_interactions = 0

        # tabular file reporting the interactions in the protein-protein interaction database
        self.input_file_location = self.configuration.get("FILE")

        # map associating database-specific identifiers with UniProt IDs
        self.id_map = {}

        # Parse configuration specifying the specific column designations of relevant entries
        # Default to common designations
        self.col_interactor_a = self.configuration.get(
            "COLUMNS", {}).get("ID INTERACTOR A")
        self.col_interactor_b = self.configuration.get(
            "COLUMNS", {}).get("ID INTERACTOR B")
        self.col_alt_interactor_a = self.configuration.get(
            "COLUMNS", {}).get("ALTERNATIVE ID INTERACTOR A")
        self.col_alt_interactor_b = self.configuration.get(
            "COLUMNS", {}).get("ALTERNATIVE ID INTERACTOR B")
        self.col_tax_id_a = self.configuration.get("COLUMNS",
                                                   {}).get("TAXONOMY ID INTERACTOR A")
        self.col_tax_id_b = self.configuration.get("COLUMNS",
                                                   {}).get("TAXONOMY ID INTERACTOR B")
        self.col_detection_method = self.configuration.get(
            "COLUMNS", {}).get("DETECTION METHOD")
        self.col_interaction_type = self.configuration.get(
            "COLUMNS", {}).get("INTERACTION TYPE")

        # Parse configuration specifying the specific column designations of confidence scores
        # if a confidence core threshold is applied
        if self.configuration.get("SCORE") is not None:
            self.col_score = self.configuration.get("COLUMNS", {}).get(
                "SCORE", "Confidence value(s)")
        else:
            self.col_score = None

        # Parse configuration specifying the delimiter of cells
        self.cell_separator = self.configuration.get("CELL DELIMITER", "\t")

        # Parse configuration specifying the delimiter of multiple entries within cells
        self.entry_separator = self.configuration.get("ENTRY DELIMITER", "|")

        # Parse configuration specifying a prefix to be excluded from the string representation
        # of a score's numerical value
        self.score_prefix = self.configuration.get("SCORE PREFIX")

        # Parse configuration specifying a factor to be removed from a score,
        # e.g. STRING scores are multiplied by 1000
        self.score_scale = float(self.configuration.get("SCORE SCALE", "1.0"))

        # Parse configuration specifying a prefix designating a required type of ID reported
        # If no file to derive a mapping between database-specific identifiers and UniProt IDs from is specified
        # default to "uniprotkb", otherwise the lower case database name, e.g. BioGRID to biogrid
        if self.configuration.get("ID MAP"):
            self.id_prefix = self.configuration.get("ID PREFIX")
        else:
            self.id_prefix = self.configuration.get("ID PREFIX", "uniprotkb")

        # Use a file specifying conversions between database-specific and UniProt IDs
        if self.configuration.get("ID MAP"):
            # download the file from a specified URL if it is not available
            map_available = os.path.exists(self.configuration["ID MAP"]["FILE"])
            if not map_available and self.configuration["ID MAP"].get("URL"):
                map_available = self.download(
                    self.configuration["ID MAP"]["URL"],
                    self.configuration["ID MAP"]["FILE"])
            if map_available:
                self.compute_map(configuration["ID MAP"]["FILE"],
                                 self.configuration["ID MAP"].get("ID", self.name))
        else:
            map_available = True

        # download the tabular file reporting interactions in the database if it is not available
        data_available = os.path.exists(self.configuration.get("FILE"))
        if not data_available and self.configuration.get("URL"):
            data_available = self.download(self.configuration["URL"],
                                           self.configuration["FILE"])

        # retrieve protein-protein interactions from the file
        if map_available and data_available:
            self.parse()

        # if a log file exists,
        # report how many interactions were obtained and discarded
        if self.logger:
            num_all_interactions = sum([
                self.num_interactions, self.discarded_score,
                self.discarded_unmapped, self.discarded_methods_types
            ])

            num_unmapped_proteins = len(self.unmapped_proteins)
            if self.interactions:
                num_mapped_proteins = len(
                    set.union(*[
                        set(interactants) for interactants in self.interactions
                    ]))
                num_proteins = num_mapped_proteins + num_unmapped_proteins

                # report to a log file how many interactions were discarded
                # because they did not meet requirements concerning detection methods or interaction types
                self.logger.info(
                    "{}: INTERACTIONS DISCARDED: DETECTION METHODS, INTERACTION TYPES: {}/{} ({}%)"
                    .format(
                        self.name, self.discarded_methods_types,
                        num_all_interactions, self.discarded_methods_types /
                        num_all_interactions * 100))

                # report to a log file how many interactions were discarded
                # because UniProt IDs could be not be determined for both interacting proteins
                self.logger.info(
                    "{}: INTERACTIONS DISCARDED: UNMAPPED PROTEIN IDENTIFIERS: {}/{} ({}%)"
                    .format(
                        self.name, self.discarded_unmapped,
                        num_all_interactions,
                        self.discarded_unmapped / num_all_interactions * 100))

                # report to a log file how many interactions were discarded
                # because their score was below the threshold
                self.logger.info(
                    "{}: INTERACTIONS DISCARDED: SCORE: {}/{} ({}%)".format(
                        self.name, self.discarded_score, num_all_interactions,
                        self.discarded_score / num_all_interactions * 100))

                # report to a log file how many proteins of those represented among the interactions
                # in the file could not be mapped be mapped to UniProt IDs,
                # assuming that there are no synonyms among these proteins
                self.logger.info(
                    "{}: PROTEINS UNMAPPED TO UNIPROT IDENTIFIERS: {}/{} ({}%)".
                    format(self.name, num_unmapped_proteins, num_proteins,
                           num_unmapped_proteins / num_proteins * 100))
            else:
                # report to a log file that no interactions were retrieved
                self.logger.info("{}:\tNO INTERACTION DATA AVAILABLE".format(
                    self.name))
                progressbar.ProgressBar(
                    max_value=progressbar.UnknownLength).finish()

    def verify_interaction_entry(self, row):
        """
        Given a row of the tabular file, determine if the corresponding cells contain 
        required PSI-MI IDs for the detection method and interaction type of an interaction
        """
        if self.configuration.get("DETECTION METHOD") and self.configuration.get("TYPE"):
            return (any((get_psimi(term) in self.configuration["DETECTION METHOD"]
                    for term in row[self.col_detection_method].split(self.cell_separator))) and
                any((get_psimi(term) in self.configuration["TYPE"]
                    for term in row[self.col_interaction_type].split(self.cell_separator))))
        elif self.configuration.get("DETECTION METHOD"):
            return any((get_psimi(term) in self.configuration["DETECTION METHOD"]
                    for term in row[self.col_detection_method].split(self.cell_separator)))
        elif self.configuration.get("TYPE"):
            return any((get_psimi(term) in self.configuration["TYPE"]
                    for term in row[self.col_interaction_type].split(self.cell_separator)))
        else:
            return True

    def verify_taxid(self, row):
        """
        Given a row of the tabular file, determine if the corresponding cells for both interacting proteins 
        contain taxonomic identifiers for human (9606) proteins 
        """
        if not self.col_tax_id_a and not self.col_tax_id_b:
            return True
        # "-" signifies an empty cell
        if row[self.col_tax_id_a].strip() == "-" or row[
                self.col_tax_id_b].strip() == "-":
            return False
        return (any(
            (int(entry.split(":")[1].split("(")[0]) == 9606
             for entry in row[self.col_tax_id_a].split(self.cell_separator))
        ) and any(
            (int(entry.split(":")[1].split("(")[0]) == 9606
             for entry in row[self.col_tax_id_b].split(self.cell_separator))))

    def compute_map(self, map_file_location, map_id):
        """Assemble a dictionary associating database-specific with UniProt IDs"""
        # obtain the identifier of the database in the file specifying the ID mapping
        # from the configuration
        self.id_map = {}
        with open(map_file_location) as map_file:
            reader = csv.reader(map_file, delimiter="\t")
            for line in reader:
                if line[1] == map_id:
                    self.id_map[line[2]] = line[0]

    def parse(self):
        """Obtain protein-protein interactions from the respective file"""
        print("{}:".format(self.name))
        with open(self.input_file_location) as input_file:
            # specify which columns of the file to import into a pandas DataFrame
            # column designations set appropriately during class initialization
            cols_to_use = [self.col_interactor_a, self.col_interactor_b]
            for col in [
                    self.col_alt_interactor_a, self.col_alt_interactor_b,
                    self.col_detection_method, self.col_interaction_type,
                    self.col_tax_id_a, self.col_tax_id_b, self.col_score
            ]:
                if col:
                    cols_to_use.append(col)

                # include columns for additional confidence values (STRING)
                for add_score in self.configuration.get("SCORES", {}):
                    cols_to_use.append(add_score["COLUMN"])

            file_data_df = pd.read_csv(
                input_file,
                # a numeric column designation applies if there is no header (MINT)
                header=None if type(self.col_interactor_a) == int else 0,
                sep=self.cell_separator,
                # simplify handling of specific values by asserting that they are strings
                dtype=str,
                low_memory=True,
                usecols=cols_to_use)
            bar = progressbar.ProgressBar(max_value=file_data_df.shape[0])
            bar.update(0)
            for i, row in file_data_df.iterrows():
                # verify that interaction is reported to occur between human proteins,
                # as opposed to for example, interactions between human and bacterial proteins
                # proceed with next interaction if this is not the case
                if not self.verify_taxid(row):
                    bar.update(i + 1)
                    continue

                # verify that the interaction is annotated with required identifiers
                # for detection methods and interaction types if applicable
                # proceed with next interaction if this is not the case
                if not self.verify_interaction_entry(row):
                    self.discarded_methods_types += 1
                    bar.update(i + 1)
                    continue
                
                # if confidence scores for interactions are considered,
                # include this interactions' confidence score in the distribution
                if self.col_score:
                    for score in row[self.col_score].split(
                            self.entry_separator):
                        # if applicable, remove the score prefix
                        if self.score_prefix:
                            if not score.split(":")[0] == self.score_prefix:
                                continue
                            score = float(
                                score.split(":")[1].replace("(inferred)", ""))
                        else:
                            score = float(score)
                        # account for the a multiplied score
                        score /= self.score_scale
                        self.scores.add(score)
                        break
                    else:
                        score = None
                else:
                    score = None

                # retrieve UniProt IDs for both interacting proteins
                # considering alternative identifiers
                ids = [None, None]
                for j, column in enumerate(
                    ((self.col_interactor_a, self.col_alt_interactor_a),
                     (self.col_interactor_b, self.col_alt_interactor_b))):
                    # search primary IDs
                    if self.id_prefix:
                        if row[column[0]].split(":")[0] == self.id_prefix:
                            ids[j] = row[column[0]].split(":")[1]
                    else:
                        ids[j] = row[column[0]]

                    if self.id_map:
                        ids[j] = self.id_map.get(ids[j]) 

                    # search alternate IDs, if available
                    if ids[j] is None and column[1]:
                        for identifier in row[column[1]].split(
                                self.entry_separator):
                            if self.id_prefix:
                                if identifier.split(":")[0] == self.id_prefix:
                                    ids[j] = identifier.split(":")[1]
                            else:
                                ids[j] = identifier

                            if self.id_map:
                                ids[j] = self.id_map.get(ids[j])
                                    
                            if ids[j] is not None:
                                break
      
                    # save primary identifier of proteins not successfully mapped to UniProt IDs
                    # to report amount to log file subsequently
                    if ids[j] is None:
                        self.unmapped_proteins.add(row[column[0]])

                # discard interaction if at least one interacting protein was not successfully mapped
                # to its respective UniProt ID
                if None in ids:
                    self.discarded_unmapped += 1
                    bar.update(i + 1)
                    continue

                # exclude self interactions
                if ids[0] == ids[1]:
                    bar.update(i + 1)
                    continue
                
                # discard interactions not meeting the (primary) score threshold
                if score and score < self.configuration.get("SCORE", 0.0):
                    self.discarded_score += 1
                    bar.update(i + 1)
                    continue

                # discard interactions not meeting additional score thresholds
                if not all((float(row[score["COLUMN"]]) >= score["THRESHOLD"]
                        for score in self.configuration.get("SCORES", {}))):
                    self.discarded_score += 1
                    bar.update(i + 1)
                    continue

                # use the largest confidence score reported for any interaction
                # between both interacting proteins
                # effectively, disregard directionality
                interactants = frozenset([ids[0], ids[1]])
                if interactants in self.interactions and self.interactions[
                        interactants] is not None:
                    score = max(score, self.interactions[interactants])

                self.interactions[interactants] = score
                self.num_interactions += 1
                bar.update(i + 1)
            bar.finish()

    def plot_score_distribution(self, file_name):
        """
        Employ the specialized class "Score" 
        to export a histogram of the confidence score distribution
        to file_name
        """
        if self.scores and file_name:
            self.scores.plot(self.name, file_name,
                             self.configuration.get("SCORE", 0.0))

    def download(self, url, file_name, chunk_size=1048576):
        """
        Download a file from a specified URL and save it at file_name
        Decompress the file if the URL implies compression (zip, gzip)
        Return True if the file is successfully provided at file_name, False otherwise
        """
        # create the directory implied by file_name to store the downloaded file
        if not os.path.exists(os.path.dirname(file_name)):
            pathlib.Path(os.path.dirname(file_name)).mkdir(parents=True,
                                                           exist_ok=True)

        # if file_name implies that the downloaded file is compressed,
        # define a (temporary) file location to download the compressed file to and extract it from
        if url.endswith(".zip") or url.endswith(".gz"):
            file_loc = os.path.join(os.path.dirname(file_name),
                                    os.path.basename(url))

        if url.startswith("http://") or url.startswith("https://"):
            # Request the file from a HTTP server
            response = requests.get(url)
            # verify that the file could be accessed successfully
            if response.status_code == 200:
                if self.logger:
                    self.logger.info("{}: DOWNLOADED {} FROM {}".format(
                        self.name, file_name, url))
            else:
                if self.logger:
                    self.logger.critical(
                        "{}:\tDOWNLOAD OF {} FROM {} FAILED: {}".format(
                            self.name, file_name, url,
                            response.status_code))
                return False

            # download the file
            # to a temporary file location, if its compressed
            # and to the final file location (file_name) otherwise
            if url.endswith(".zip") or url.endswith(".gz"):
                file_loc = os.path.join(os.path.dirname(file_name),
                                        os.path.basename(url))
                with open(file_loc, "wb") as downloaded_file:
                    for chunk in response.iter_content(chunk_size=chunk_size):
                        downloaded_file.write(chunk)
            else:
                with open(file_name, "wb") as downloaded_file:
                    for chunk in response.iter_content(chunk_size=chunk_size):
                        downloaded_file.write(chunk)

        elif url.startswith("ftp://"):
            # Obtain the file from an FTP server
            ftp = ftplib.FTP(urlparse(url).netloc)
            # Login to the server
            response_message = ftp.login()
            path = urlparse(url).path
            # change to the directory containing the file
            ftp.cwd(os.path.dirname(path))

            # verify that the file could be accessed successfully
            if response_message.startswith("2"):
                if self.logger:
                    self.logger.info("{}: DOWNLOADED {} FROM {}".format(
                        self.name, file_name, url))
            else:
                if self.logger:
                    self.logger.critical(
                        "{}:\tDOWNLOAD OF {} FROM {} FAILED: {}".format(
                            self.name, file_name, url,
                            response_message))
                return False

            # download the file
            # to a temporary file location, if its compressed
            # and to the final file location (file_name) otherwise
            if url.endswith(".zip") or url.endswith(".gz"):
                with open(file_loc, "wb") as downloaded_file:
                    ftp.retrbinary("RETR {}".format(os.path.basename(path)),
                                   downloaded_file.write,
                                   blocksize=chunk_size)
            else:
                with open(file_name, "wb") as downloaded_file:
                    ftp.retrbinary("RETR {}".format(os.path.basename(path)),
                                   downloaded_file.write,
                                   blocksize=chunk_size)
            ftp.quit()

        if url.endswith(".zip") or url.endswith(".gz"):
            # decompress the file
            if url.endswith(".zip"):
                try:
                    with zipfile.ZipFile(file_loc) as decompressed_archive:
                        decompressed_archive.extract(
                            os.path.basename(file_name),
                            os.path.dirname(file_name))
                except:
                    return False
            else:
                try:
                    with open(file_name, "wb") as out_file:
                        with gzip.open(file_loc, "rb") as decompressed_file:
                            chunk = decompressed_file.read(chunk_size)
                            while chunk:
                                out_file.write(chunk)
                                chunk = decompressed_file.read(chunk_size)
                except:
                    return False

            # remove the compressed file
            os.remove(file_loc)
        return True

    class Scores():
        """
        Class to store confidence score distribution of a protein-protein interaction database memory efficiently.
        Given a fixed number of equal length intervalls within a range, individual scores are aggregated in intervalls.
        Intervalls are represented by counters representing the number of scores within this intervall.
        """

        def __init__(self, data_range=(0.0, 1.0), values=100):
            self.num_bins = values
            self.bins = [0 for i in range(self.num_bins + 1)]
            self.range = data_range[1] - data_range[0]
            self.min = data_range[0]
            self.max = data_range[1]
            self.n = 0
            self.bin_width = self.range / self.num_bins

        def __bool__(self):
            """Return true if distribution contains observations"""
            return any(self.bins)

        def add(self, value):
            """Add individual score to distribution"""
            if value is not None:
                self.bins[int((value - self.min) / self.bin_width +
                              0.5 * self.bin_width)] += 1
                self.n += 1

        def add_from(self, values):
            """Add multiple individual scores to distribution"""
            for value in values:
                self.add(value)

        def plot(self, name, file_name, score_threshold, num_plot_bins=20):
            """
            Export a histogram of the represented distribution to file_name.
            """
            plot_bin_width = self.range / num_plot_bins
            plt.xlabel(r"{} score".format(name))
            plt.ylabel(r"interactions [%]")
            plt.xticks(np.arange(0.0, 1.1, 0.1))
            plt.yticks([
                self.n * s
                for s in np.arange(self.min, self.max +
                                   self.range * 0.2, self.range * 0.2)
            ])
            plt.gca().yaxis.set_major_formatter(
                ticker.PercentFormatter(decimals=0,
                                        symbol=None,
                                        xmax=self.n,
                                        is_latex=True))
            plt.xlim(self.min - 0.5 * plot_bin_width,
                     self.max + 0.5 * plot_bin_width)
            plt.ylim(0, self.n)
            plt.hist(np.arange(self.min, self.max + self.bin_width,
                               self.bin_width),
                     weights=self.bins,
                     bins=np.arange(self.min + plot_bin_width, self.max +
                                    2 * plot_bin_width, plot_bin_width) - 0.5 *
                     (plot_bin_width),
                     edgecolor="k",
                     linewidth=1.0,
                     color="w")

            if score_threshold:
                # plot a vertical line indicating the confidence score threshold applied to the data
                plt.gca().axvline(score_threshold,
                                  color="k",
                                  linestyle="dashed",
                                  linewidth=1.0)
            plt.savefig("{0}.{2}{1}".format(*os.path.splitext(file_name),
                                            name.lower()),
                        dpi=DPI)
            plt.close()
