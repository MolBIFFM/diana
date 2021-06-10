import bisect
import itertools
import os
import math
import statistics
import concurrent.futures

import networkx as nx
import scipy.stats
import pandas as pd

from pipeline.configuration import data
from pipeline.utilities import fetch, mitab, modularization, correction


class ProteinInteractionNetwork(nx.Graph):
    def __init__(self):
        super().__init__()

    def annotate_proteins(self):
        accessions, entry_gene_name, entry_protein_name = [], {}, {}
        rec_name = False
        for line in fetch.txt(data.UNIPROT_SWISSPROT):
            if not line.strip():
                continue
            if line.split(maxsplit=1)[0] == "AC":
                if len(line.split(maxsplit=1)) == 1:
                    continue

                accessions.extend(line.split(maxsplit=1)[1].rstrip(";").split("; "))

            elif line.split(maxsplit=1)[0] == "GN":
                if len(line.split(maxsplit=1)) == 1:
                    continue

                for entry in line.split(maxsplit=1)[1].rstrip(";").split("; "):
                    if "=" in entry:
                        entry_gene_name[entry.split("=")[0]] = (
                            entry.split("=")[1].split("{")[0].rstrip()
                        )

            elif line.split(maxsplit=1)[0] == "DE":
                if len(line.split(maxsplit=1)) == 1:
                    continue

                if line.split(maxsplit=1)[1].split(":", 1)[0] == "RecName":
                    entries = (
                        line.split(maxsplit=1)[1]
                        .split(":", 1)[1]
                        .lstrip()
                        .rstrip(";")
                        .split("; ")
                    )
                    rec_name = True
                elif line.split(maxsplit=1)[1].split(":", 1)[0] == "AltName":
                    entries = []
                    rec_name = False
                elif line.split(maxsplit=1)[1].split(":", 1)[0] == "Contains":
                    entries = []
                    rec_name = False
                elif line.split(maxsplit=1)[1].split(":", 1)[0] == "Flags":
                    entries = []
                    rec_name = False
                elif rec_name:
                    entries = line.split(maxsplit=1)[1].rstrip(";").split("; ")

                for entry in entries:
                    if "=" in entry:
                        if entry.split("=")[0] not in entry_protein_name:
                            entry_protein_name[entry.split("=")[0]] = (
                                entry.split("=")[1].split("{")[0].rstrip()
                            )

            elif line == "//":
                for protein in tuple(self):
                    if protein.split("-")[0] in accessions:
                        if "-" in protein and not protein.split("-")[1].isnumeric():
                            nx.relabel_nodes(
                                self, {protein: protein.split("-")[0]}, copy=False
                            )
                            protein = protein.split("-")[0]

                        if accessions.index(protein.split("-")[0]) == 0:
                            self.nodes[protein]["gene"] = entry_gene_name.get(
                                "Name", "NA"
                            )
                            self.nodes[protein]["protein"] = entry_protein_name.get(
                                "Full", "NA"
                            )
                        else:
                            nx.relabel_nodes(self, {protein: accessions[0]}, copy=False)
                            self.nodes[accessions[0]]["gene"] = entry_gene_name.get(
                                "Name", "NA"
                            )
                            self.nodes[accessions[0]][
                                "protein"
                            ] = entry_protein_name.get("Full", "NA")

                accessions.clear()
                entry_gene_name.clear()
                entry_protein_name.clear()

                rec_name = False

    def remove_unannotated_proteins(self):
        self.remove_nodes_from(
            [node for node in self if not self.nodes[node].get("protein")]
        )

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

            genes.update(
                [
                    str(gene_accession)
                    for gene_accession in gene_accession_format(
                        row[gene_accession_column]
                    )
                ]
            )

        accessions, entry_gene_name, entry_protein_name = [], {}, {}
        rec_name, taxon_id = False, 0
        for line in fetch.txt(data.UNIPROT_SWISSPROT):
            if not line.strip():
                continue
            if line.split(maxsplit=1)[0] == "AC":
                if len(line.split(maxsplit=1)) == 1:
                    continue

                accessions.extend(line.split(maxsplit=1)[1].rstrip(";").split("; "))

            elif line.split(maxsplit=1)[0] == "GN":
                if len(line.split(maxsplit=1)) == 1:
                    continue

                for entry in line.split(maxsplit=1)[1].rstrip(";").split("; "):
                    if "=" in entry:
                        entry_gene_name[entry.split("=")[0]] = (
                            entry.split("=")[1].split("{")[0].rstrip()
                        )

            elif line.split(maxsplit=1)[0] == "DE":
                if len(line.split(maxsplit=1)) == 1:
                    continue

                if line.split(maxsplit=1)[1].split(":", 1)[0] == "RecName":
                    entries = (
                        line.split(maxsplit=1)[1]
                        .split(":", 1)[1]
                        .lstrip()
                        .rstrip(";")
                        .split("; ")
                    )
                    rec_name = True
                elif line.split(maxsplit=1)[1].split(":", 1)[0] == "AltName":
                    entries = []
                    rec_name = False
                elif line.split(maxsplit=1)[1].split(":", 1)[0] == "Contains":
                    entries = []
                    rec_name = False
                elif line.split(maxsplit=1)[1].split(":", 1)[0] == "Flags":
                    entries = []
                    rec_name = False
                elif rec_name:
                    entries = line.split(maxsplit=1)[1].rstrip(";").split("; ")

                for entry in entries:
                    if "=" in entry:
                        if entry.split("=")[0] not in entry_protein_name:
                            entry_protein_name[entry.split("=")[0]] = (
                                entry.split("=")[1].split("{")[0].rstrip()
                            )

            elif line.split(maxsplit=1)[0] == "OX":
                if len(line.split(maxsplit=1)) == 1:
                    continue

                if (
                    line.split(maxsplit=1)[1].split(";")[0].split("=")[0]
                    == "NCBI_TaxID"
                ):
                    if (
                        line.split(maxsplit=1)[1]
                        .split(";")[0]
                        .split("=")[1]
                        .split("{")[0]
                        .isnumeric()
                    ):
                        taxon_id = int(
                            line.split(maxsplit=1)[1]
                            .split(";")[0]
                            .split("=")[1]
                            .split("{")[0]
                        )

            elif line == "//":
                if (
                    taxon_id == taxon_identifier
                    and entry_gene_name.get("Name") in genes
                ):
                    self.add_node(accessions[0])
                    self.nodes[accessions[0]]["gene"] = entry_gene_name.get(
                        "Name", "NA"
                    )
                    self.nodes[accessions[0]]["protein"] = entry_protein_name.get(
                        "Full", "NA"
                    )

                accessions.clear()
                entry_gene_name.clear()
                entry_protein_name.clear()

                rec_name, taxon_id = False, 0

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
                usecols=[protein_accession_column]
                + [column for column in (position_column,) if column]
                + replicates,
                dtype={
                    protein_accession_column: str,
                    **{column: str for column in (position_column,) if column},
                    **{column: float for column in replicates},
                },
            )
        else:
            table = pd.read_table(
                file_name,
                header=header,
                usecols=[protein_accession_column]
                + [column for column in (position_column,) if column]
                + replicates,
                dtype={
                    protein_accession_column: str,
                    **{column: str for column in (position_column,) if column},
                    **{column: float for column in replicates},
                },
            )

        proteins = {}
        for _, row in table.iterrows():
            if pd.isna(row[protein_accession_column]):
                continue

            if position_column and pd.isna(row[position_column]):
                continue

            if replicates:
                measurements = [
                    row[repl] for repl in replicates if not pd.isna(row[repl])
                ]

                if len(measurements) >= min(num_replicates, len(replicates)):
                    protein_accessions = [
                        str(protein_accession)
                        for protein_accession in protein_accession_format(
                            row[protein_accession_column]
                        )
                    ]

                    if position_column:
                        positions = [
                            int(position)
                            for position in position_format(row[position_column])
                        ]
                    else:
                        positions = []

                    if len(protein_accessions) != len(positions):
                        if len(protein_accessions) > len(positions):
                            positions.extend(
                                [
                                    0
                                    for _ in range(
                                        len(positions), len(protein_accessions)
                                    )
                                ]
                            )
                        else:
                            positions = positions[: len(protein_accessions)]

                    for protein_accession, position in zip(
                        protein_accessions, positions
                    ):
                        if "-" in protein_accession:
                            protein, isoform = protein_accession.split("-")
                        else:
                            protein, isoform = protein_accession, "0"

                        if protein not in proteins:
                            proteins[protein] = {}

                        if isoform not in proteins[protein]:
                            proteins[protein][isoform] = []

                        bisect.insort(
                            proteins[protein][isoform],
                            [
                                position,
                                convert_measurement(combine_replicates(measurements)),
                            ],
                        )
            else:
                protein_accessions = [
                    str(protein_accession)
                    for protein_accession in protein_accession_format(
                        row[protein_accession_column]
                    )
                ]

                for protein_accession in protein_accessions:
                    if "-" in protein_accession:
                        protein, isoform = protein_accession.split("-")
                    else:
                        protein, isoform = protein_accession, "0"

                    if protein not in proteins:
                        proteins[protein] = {}

                    if isoform not in proteins[protein]:
                        proteins[protein][isoform] = []

        swissprot_proteins, primary_accession, gene_name, protein_name = {}, {}, {}, {}

        accessions, entry_gene_name, entry_protein_name = [], {}, {}
        rec_name = False
        for line in fetch.txt(data.UNIPROT_SWISSPROT):
            if not line.strip():
                continue

            if line.split(maxsplit=1)[0] == "AC":
                if len(line.split(maxsplit=1)) == 1:
                    continue

                accessions.extend(line.split(maxsplit=1)[1].rstrip(";").split("; "))

            elif line.split(maxsplit=1)[0] == "GN":
                if len(line.split(maxsplit=1)) == 1:
                    continue

                for entry in line.split(maxsplit=1)[1].rstrip(";").split("; "):
                    if "=" in entry:
                        entry_gene_name[entry.split("=")[0]] = (
                            entry.split("=")[1].split("{")[0].rstrip()
                        )

            elif line.split(maxsplit=1)[0] == "DE":
                if len(line.split(maxsplit=1)) == 1:
                    continue

                if line.split(maxsplit=1)[1].split(":", 1)[0] == "RecName":
                    entries = (
                        line.split(maxsplit=1)[1]
                        .split(":", 1)[1]
                        .lstrip()
                        .rstrip(";")
                        .split("; ")
                    )
                    rec_name = True
                elif line.split(maxsplit=1)[1].split(":", 1)[0] == "AltName":
                    entries = []
                    rec_name = False
                elif line.split(maxsplit=1)[1].split(":", 1)[0] == "Contains":
                    entries = []
                    rec_name = False
                elif line.split(maxsplit=1)[1].split(":", 1)[0] == "Flags":
                    entries = []
                    rec_name = False
                elif rec_name:
                    entries = line.split(maxsplit=1)[1].rstrip(";").split("; ")

                for entry in entries:
                    if "=" in entry:
                        if entry.split("=")[0] not in entry_protein_name:
                            entry_protein_name[entry.split("=")[0]] = (
                                entry.split("=")[1].split("{")[0].rstrip()
                            )

            elif line == "//":
                for i, accession in enumerate(accessions):
                    if accession in proteins:
                        swissprot_proteins[accession] = proteins[accession]
                        gene_name[accession] = entry_gene_name.get("Name", "NA")
                        protein_name[accession] = entry_protein_name.get("Full", "NA")

                        if i > 0:
                            primary_accession[accession] = accessions[0]

                accessions.clear()
                entry_gene_name.clear()
                entry_protein_name.clear()

                rec_name = False

        proteins = swissprot_proteins

        for protein in primary_accession:
            if primary_accession[protein] in proteins:
                for isoform in proteins[protein]:
                    for primary_isoform in proteins[primary_accession[protein]]:
                        for i in range(len(proteins[protein][isoform])):
                            for j in range(
                                len(
                                    proteins[primary_accession[protein]][
                                        primary_isoform
                                    ]
                                )
                            ):
                                if not proteins[protein][isoform][i]:
                                    continue

                                if (
                                    isoform == primary_isoform
                                    and proteins[protein][isoform][i][1]
                                    == proteins[primary_accession[protein]][
                                        primary_isoform
                                    ][j][1]
                                ):
                                    if not proteins[primary_accession[protein]][
                                        primary_isoform
                                    ][j][0]:
                                        proteins[primary_accession[protein]][
                                            primary_isoform
                                        ][j][0] = proteins[protein][isoform][i][0]

                                    proteins[protein][isoform][i] = []

                    if not isoform in proteins[primary_accession[protein]]:
                        proteins[primary_accession[protein]][isoform] = []

                    for measurement in proteins[protein][isoform]:
                        if measurement:
                            bisect.insort(
                                proteins[primary_accession[protein]][isoform],
                                measurement,
                            )
            else:
                proteins[primary_accession[protein]] = proteins[protein]
                gene_name[primary_accession[protein]] = gene_name[protein]
                protein_name[primary_accession[protein]] = protein_name[protein]

        for protein in primary_accession:
            del proteins[protein]

        for protein in proteins:
            for isoform in proteins[protein]:
                if isoform != "0":
                    self.add_node("-".join([protein, isoform]))
                    self.nodes["-".join([protein, isoform])]["gene"] = gene_name.get(
                        protein, "NA"
                    )
                    self.nodes["-".join([protein, isoform])][
                        "protein"
                    ] = protein_name.get(protein, "NA")
                else:
                    self.add_node(protein)
                    self.nodes[protein]["gene"] = gene_name.get(protein, "NA")
                    self.nodes[protein]["protein"] = protein_name.get(protein, "NA")

                proteins[protein][isoform] = sorted(
                    sorted(
                        proteins[protein][isoform],
                        key=lambda item: item[1],
                        reverse=True,
                    )[:num_sites]
                )

                for i in range(len(proteins[protein][isoform])):
                    if isoform != "0":
                        self.nodes["-".join([protein, isoform])][
                            "{} {} {}".format(time, ptm, i + 1)
                        ] = proteins[protein][isoform][i][1]
                    else:
                        self.nodes[protein][
                            "{} {} {}".format(time, ptm, i + 1)
                        ] = proteins[protein][isoform][i][1]

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
                if len(change.split(" ")) == 3 and change.split(" ")[0] == str(time)
            )
        )

    def get_sites(self, time, ptm):
        return max(
            int(change.split(" ")[2])
            for protein in self
            for change in self.nodes[protein]
            if len(change.split(" ")) == 3
            and change.split(" ")[0] == str(time)
            and change.split(" ")[1] == ptm
        )

    def set_post_translational_modification(self):
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

    def get_changes(self, time, ptm, combine_sites=statistics.mean):
        changes = []
        for protein in self:
            sites = [
                self.nodes[protein][change]
                for change in self.nodes[protein]
                if len(change.split(" ")) == 3
                and change.split(" ")[0] == str(time)
                and change.split(" ")[1] == ptm
            ]

            if sites:
                changes.append(combine_sites(sites))

        return changes

    def get_z_score_range(
        self,
        time,
        post_translational_modification,
        thresholds,
        combine_sites=statistics.mean,
    ):
        changes = self.get_changes(time, post_translational_modification, combine_sites)
        mean = statistics.mean(changes)
        stdev = statistics.stdev(changes, xbar=mean)
        return (thresholds[0] * stdev + mean, thresholds[1] * stdev + mean)

    def get_z_score(
        self,
        time,
        post_translational_modification,
        change,
        combine_sites=statistics.mean,
    ):
        changes = self.get_changes(time, post_translational_modification, combine_sites)
        mean = statistics.mean(changes)
        stdev = statistics.stdev(changes, xbar=mean)
        return (change - mean) / stdev

    def get_proportion_range(
        self,
        time,
        post_translational_modification,
        thresholds,
        combine_sites=statistics.mean,
    ):
        changes = sorted(
            self.get_changes(time, post_translational_modification, combine_sites)
        )
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
        post_translational_modification,
        change,
        combine_sites=statistics.mean,
    ):
        changes = self.get_changes(time, post_translational_modification, combine_sites)
        return (
            min(
                len([c for c in changes if c <= change]),
                len([c for c in changes if c >= change]),
            )
            / len(changes)
        )

    def set_changes(
        self,
        combine_sites=statistics.mean,
        changes=(-1.0, 1.0),
        get_range=lambda time, post_translational_modification, changes, combine_sites: changes,
    ):
        times = self.get_times()
        post_translational_modifications = {
            time: self.get_post_translational_modifications(time) for time in times
        }

        mid_range = {
            time: {
                post_translational_modification: get_range(
                    time,
                    post_translational_modification,
                    changes,
                    combine_sites,
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
                        and change.split(" ")[0] == str(time)
                        and change.split(" ")[1] == post_translational_modification
                    ]

                    if sites:
                        ptms[post_translational_modification] = combine_sites(sites)

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
                else:
                    self.nodes[protein]["change {}".format(time)] = "mid"

    def get_databases(self):
        return sorted(
            set(
                database for edge in self.edges for database in self.edges[edge]
            ).intersection({"BioGRID", "CORUM", "IntAct", "Reactome", "STRING"})
        )

    def set_edge_weights(
        self,
        weight=lambda confidence_scores: int(bool(confidence_scores)),
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

    def add_proteins_from_biogrid(
        self,
        experimental_system=[],
        experimental_system_type=["physical"],
        taxon_identifier=9606,
        multi_validated_physical=False,
    ):
        uniprot = {}
        for row in fetch.tabular_txt(
            data.UNIPROT_ID_MAP.format(
                organism=data.ORGANISM[taxon_identifier][data.UNIPROT_ID_MAP],
                taxon_identifier=taxon_identifier,
            ),
            delimiter="\t",
            usecols=[0, 1, 2],
        ):
            if row[1] == "BioGRID":
                uniprot[int(row[2])] = row[0]

        for row in fetch.tabular_txt(
            data.BIOGRID_ID_MAP_ZIP_ARCHIVE,
            file=data.BIOGRID_ID_MAP,
            delimiter="\t",
            header=20,
            usecols=[
                "BIOGRID_ID",
                "IDENTIFIER_VALUE",
                "IDENTIFIER_TYPE",
            ],
        ):
            if row["IDENTIFIER_TYPE"] in ("UNIPROT-ACCESSION", "UNIPROT-ISOFORM"):
                uniprot[int(row["BIOGRID_ID"])] = row["IDENTIFIER_VALUE"]

        for row in fetch.tabular_txt(
            data.BIOGRID_MV_PHYSICAL_ZIP_ARCHIVE
            if multi_validated_physical
            else data.BIOGRID_ZIP_ARCHIVE,
            file=data.BIOGRID_MV_PHYSICAL
            if multi_validated_physical
            else data.BIOGRID.format(
                organism=data.ORGANISM[taxon_identifier][data.BIOGRID]
            ),
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
                row["BioGRID ID Interactor A"] in uniprot
                and row["BioGRID ID Interactor B"] in uniprot
                and row["Organism ID Interactor A"] == row["Organism ID Interactor B"]
                and (
                    not experimental_system
                    or row["Experimental System"] in experimental_system
                )
                and (
                    not experimental_system_type
                    or row["Experimental System Type"] in experimental_system_type
                )
            ):
                if (
                    uniprot[row["BioGRID ID Interactor A"]] in self
                    and uniprot[row["BioGRID ID Interactor B"]] not in self
                    and self.nodes[uniprot[row["BioGRID ID Interactor A"]]].get(
                        "protein"
                    )
                ):
                    self.add_node(uniprot[row["BioGRID ID Interactor B"]])

                elif (
                    uniprot[row["BioGRID ID Interactor A"]] not in self
                    and uniprot[row["BioGRID ID Interactor B"]] in self
                    and self.nodes[uniprot[row["BioGRID ID Interactor B"]]].get(
                        "protein"
                    )
                ):
                    self.add_node(uniprot[row["BioGRID ID Interactor A"]])

    def add_interactions_from_biogrid(
        self,
        experimental_system=[],
        experimental_system_type=["physical"],
        taxon_identifier=9606,
        multi_validated_physical=False,
    ):
        uniprot = {}
        for row in fetch.tabular_txt(
            data.UNIPROT_ID_MAP.format(
                organism=data.ORGANISM[taxon_identifier][data.UNIPROT_ID_MAP],
                taxon_identifier=taxon_identifier,
            ),
            delimiter="\t",
            usecols=[0, 1, 2],
        ):
            if row[1] == "BioGRID" and row[0] in self:
                uniprot[int(row[2])] = row[0]

        for row in fetch.tabular_txt(
            data.BIOGRID_ID_MAP_ZIP_ARCHIVE,
            file=data.BIOGRID_ID_MAP,
            delimiter="\t",
            header=20,
            usecols=[
                "BIOGRID_ID",
                "IDENTIFIER_VALUE",
                "IDENTIFIER_TYPE",
            ],
        ):
            if (
                row["IDENTIFIER_TYPE"] in ("UNIPROT-ACCESSION", "UNIPROT-ISOFORM")
                and row["IDENTIFIER_VALUE"] in self
            ):
                uniprot[int(row["BIOGRID_ID"])] = row["IDENTIFIER_VALUE"]

        for row in fetch.tabular_txt(
            data.BIOGRID_MV_PHYSICAL_ZIP_ARCHIVE
            if multi_validated_physical
            else data.BIOGRID_ZIP_ARCHIVE,
            file=data.BIOGRID_MV_PHYSICAL
            if multi_validated_physical
            else data.BIOGRID.format(
                organism=data.ORGANISM[taxon_identifier][data.BIOGRID]
            ),
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
                row["BioGRID ID Interactor A"] in uniprot
                and row["BioGRID ID Interactor B"] in uniprot
                and uniprot[row["BioGRID ID Interactor A"]]
                != uniprot[row["BioGRID ID Interactor B"]]
                and (
                    not experimental_system
                    or row["Experimental System"] in experimental_system
                )
                and (
                    not experimental_system_type
                    or row["Experimental System Type"] in experimental_system_type
                )
            ):
                self.add_edge(
                    uniprot[row["BioGRID ID Interactor A"]],
                    uniprot[row["BioGRID ID Interactor B"]],
                )
                self.edges[
                    uniprot[row["BioGRID ID Interactor A"]],
                    uniprot[row["BioGRID ID Interactor B"]],
                ]["BioGRID"] = 0.5

    def add_proteins_from_corum(self, protein_complex_purification_method=[]):
        for row in fetch.tabular_txt(
            data.CORUM_ZIP_ARCHIVE,
            file=data.CORUM,
            delimiter="\t",
            header=0,
            usecols=[
                "subunits(UniProt IDs)",
                "Protein complex purification method",
            ],
        ):
            if not protein_complex_purification_method or any(
                method in protein_complex_purification_method
                for method in [
                    entry.split("-")[1].lstrip()
                    for entry in row["Protein complex purification method"].split(";")
                ]
            ):
                for interactor_a, interactor_b in itertools.combinations(
                    row["subunits(UniProt IDs)"].split(";"), 2
                ):
                    if (
                        interactor_a in self
                        and interactor_b not in self
                        and self.nodes[interactor_a].get("protein")
                    ):
                        self.add_node(interactor_b)

                    elif (
                        interactor_a not in self
                        and interactor_b in self
                        and self.nodes[interactor_b].get("protein")
                    ):
                        self.add_node(interactor_a)

    def add_interactions_from_corum(self, protein_complex_purification_method=[]):
        for row in fetch.tabular_txt(
            data.CORUM_ZIP_ARCHIVE,
            file=data.CORUM,
            delimiter="\t",
            header=0,
            usecols=[
                "subunits(UniProt IDs)",
                "Protein complex purification method",
            ],
        ):
            if not protein_complex_purification_method or any(
                method in protein_complex_purification_method
                for method in [
                    entry.split("-")[1].lstrip()
                    for entry in row["Protein complex purification method"].split(";")
                ]
            ):
                for interactor_a, interactor_b in itertools.combinations(
                    row["subunits(UniProt IDs)"].split(";"), 2
                ):
                    if interactor_a in self and interactor_b in self:
                        self.add_edge(
                            interactor_a,
                            interactor_b,
                        )
                        self.edges[
                            interactor_a,
                            interactor_b,
                        ]["CORUM"] = 0.5

    def add_proteins_from_intact(
        self,
        interaction_detection_methods=[],
        interaction_types=[],
        mi_score=0.0,
    ):
        for row in fetch.tabular_txt(
            data.INTACT_ZIP_ARCHIVE,
            file=data.INTACT,
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

            if (
                (
                    interactor_a := mitab.get_identifier_from_namespace(
                        row["#ID(s) interactor A"], "uniprotkb"
                    )
                )
                and (
                    interactor_b := mitab.get_identifier_from_namespace(
                        row["ID(s) interactor B"], "uniprotkb"
                    )
                )
                and mitab.get_identifier_from_namespace(
                    row["Taxid interactor A"], "taxid"
                )
                == mitab.get_identifier_from_namespace(
                    row["Taxid interactor B"], "taxid"
                )
                and (
                    not interaction_detection_methods
                    or mitab.namespace_has_any_term_from(
                        row["Interaction detection method(s)"],
                        "psi-mi",
                        interaction_detection_methods,
                    )
                )
                and (
                    not interaction_types
                    or mitab.namespace_has_any_term_from(
                        row["Interaction type(s)"], "psi-mi", interaction_types
                    )
                )
                and (
                    score := mitab.get_identifier_from_namespace(
                        row["Confidence value(s)"], "intact-miscore"
                    )
                )
                and float(score) >= mi_score
            ):
                if (
                    interactor_a in self
                    and interactor_b not in self
                    and self.nodes[interactor_a].get("protein")
                ):
                    self.add_node(interactor_b)

                elif (
                    interactor_a not in self
                    and interactor_b in self
                    and self.nodes[interactor_b].get("protein")
                ):
                    self.add_node(interactor_a)

    def add_interactions_from_intact(
        self,
        interaction_detection_methods=[],
        interaction_types=[],
        mi_score=0.0,
    ):

        for row in fetch.tabular_txt(
            data.INTACT_ZIP_ARCHIVE,
            file=data.INTACT,
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

            if (
                (
                    interactor_a := mitab.get_identifier_from_namespace(
                        row["#ID(s) interactor A"], "uniprotkb"
                    )
                )
                and (
                    interactor_b := mitab.get_identifier_from_namespace(
                        row["ID(s) interactor B"], "uniprotkb"
                    )
                )
                and interactor_a != interactor_b
                and (
                    not interaction_detection_methods
                    or mitab.namespace_has_any_term_from(
                        row["Interaction detection method(s)"],
                        "psi-mi",
                        interaction_detection_methods,
                    )
                )
                and (
                    not interaction_types
                    or mitab.namespace_has_any_term_from(
                        row["Interaction type(s)"], "psi-mi", interaction_types
                    )
                )
                and (
                    score := mitab.get_identifier_from_namespace(
                        row["Confidence value(s)"], "intact-miscore"
                    )
                )
                and interactor_a in self
                and interactor_b in self
            ):
                score = float(score)

                if score >= mi_score:
                    if self.has_edge(interactor_a, interactor_b):
                        self.edges[interactor_a, interactor_b]["IntAct"] = max(
                            score,
                            self.edges[interactor_a, interactor_b].get("IntAct", 0.0),
                        )
                    else:
                        self.add_edge(interactor_a, interactor_b)
                        self.edges[interactor_a, interactor_b]["IntAct"] = score

    def add_proteins_from_reactome(
        self,
        interaction_type=[],
        interaction_context=[],
        taxon_identifier=9606,
    ):
        for row in fetch.tabular_txt(
            data.REACTOME.format(
                organism=data.ORGANISM[taxon_identifier][data.REACTOME]
            ),
            delimiter="\t",
            header=0,
            usecols=[
                "# Interactor 1 uniprot id",
                "Interactor 2 uniprot id",
                "Interaction type",
                "Interaction context",
            ],
        ):
            if (
                row["# Interactor 1 uniprot id"].split(":")[0] == "uniprotkb"
                and row["Interactor 2 uniprot id"].split(":")[0] == "uniprotkb"
            ):
                interactor_a = row["# Interactor 1 uniprot id"].split(":")[1]
                interactor_b = row["Interactor 2 uniprot id"].split(":")[1]

                if (
                    not interaction_type or row["Interaction type"] in interaction_type
                ) and (
                    not interaction_context
                    or row["Interaction context"] in interaction_context
                ):
                    if (
                        interactor_a in self
                        and interactor_b not in self
                        and self.nodes[interactor_a].get("protein")
                    ):
                        self.add_node(interactor_b)

                    elif (
                        interactor_a not in self
                        and interactor_b in self
                        and self.nodes[interactor_b].get("protein")
                    ):
                        self.add_node(interactor_a)

    def add_interactions_from_reactome(
        self,
        interaction_type=[],
        interaction_context=[],
        taxon_identifier=9606,
    ):

        for row in fetch.tabular_txt(
            data.REACTOME.format(
                organism=data.ORGANISM[taxon_identifier][data.REACTOME]
            ),
            delimiter="\t",
            header=0,
            usecols=[
                "# Interactor 1 uniprot id",
                "Interactor 2 uniprot id",
                "Interaction type",
                "Interaction context",
            ],
        ):
            if (
                row["# Interactor 1 uniprot id"].split(":")[0] == "uniprotkb"
                and row["Interactor 2 uniprot id"].split(":")[0] == "uniprotkb"
            ):
                interactor_a = row["# Interactor 1 uniprot id"].split(":")[1]
                interactor_b = row["Interactor 2 uniprot id"].split(":")[1]

                if (
                    interactor_a != interactor_b
                    and (
                        not interaction_type
                        or row["Interaction type"] in interaction_type
                    )
                    and (
                        not interaction_context
                        or row["Interaction context"] in interaction_context
                    )
                    and interactor_a in self
                    and interactor_b in self
                ):

                    self.add_edge(interactor_a, interactor_b)
                    self.edges[interactor_a, interactor_b]["Reactome"] = 0.5

    def add_proteins_from_string(
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
        taxon_identifier=9606,
        physical=False,
    ):

        uniprot = {}
        for row in fetch.tabular_txt(
            data.UNIPROT_ID_MAP.format(
                organism=data.ORGANISM[taxon_identifier][data.UNIPROT_ID_MAP],
                taxon_identifier=taxon_identifier,
            ),
            delimiter="\t",
            usecols=[0, 1, 2],
        ):
            if row[1] == "STRING":
                uniprot[row[2]] = row[0]

        for row in fetch.tabular_txt(
            data.STRING_ID_MAP.format(taxon_identifier=taxon_identifier),
            delimiter="\t",
            skiprows=1,
            usecols=[0, 1, 2],
        ):
            if "BLAST_UniProt_AC" in row[2].split():
                uniprot[row[0]] = row[1]

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

        for row in fetch.tabular_txt(
            data.STRING_PHYSICAL.format(taxon_identifier=taxon_identifier)
            if physical
            else data.STRING.format(taxon_identifier=taxon_identifier),
            delimiter=" ",
            header=0,
            usecols=["protein1", "protein2"] + list(thresholds.keys()),
        ):
            if (
                row["protein1"] in uniprot
                and row["protein2"] in uniprot
                and uniprot[row["protein1"]] != uniprot[row["protein2"]]
                and all(
                    row[column] / 1000 >= thresholds[column] for column in thresholds
                )
            ):
                if (
                    uniprot[row["protein1"]] in self
                    and uniprot[row["protein2"]] not in self
                    and self.nodes[uniprot[row["protein1"]]].get("protein")
                ):
                    self.add_node(uniprot[row["protein2"]])

                elif (
                    uniprot[row["protein1"]] not in self
                    and uniprot[row["protein2"]] in self
                    and self.nodes[uniprot[row["protein2"]]].get("protein")
                ):
                    self.add_node(uniprot[row["protein1"]])

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
        taxon_identifier=9606,
        physical=False,
    ):

        uniprot = {}
        for row in fetch.tabular_txt(
            data.UNIPROT_ID_MAP.format(
                organism=data.ORGANISM[taxon_identifier][data.UNIPROT_ID_MAP],
                taxon_identifier=taxon_identifier,
            ),
            delimiter="\t",
            usecols=[0, 1, 2],
        ):
            if row[1] == "STRING" and row[0] in self:
                uniprot[row[2]] = row[0]

        for row in fetch.tabular_txt(
            data.STRING_ID_MAP.format(taxon_identifier=taxon_identifier),
            delimiter="\t",
            skiprows=1,
            usecols=[0, 1, 2],
        ):
            if "BLAST_UniProt_AC" in row[2].split() and row[1] in self:
                uniprot[row[0]] = row[1]

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

        for row in fetch.tabular_txt(
            data.STRING_PHYSICAL.format(taxon_identifier=taxon_identifier)
            if physical
            else data.STRING.format(taxon_identifier=taxon_identifier),
            delimiter=" ",
            header=0,
            usecols=["protein1", "protein2"] + list(thresholds.keys()),
        ):
            if (
                row["protein1"] in uniprot
                and row["protein2"] in uniprot
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

    def get_modules(
        self,
        module_size=35,
        combine_sizes=statistics.mean,
        weight="weight",
        algorithm=modularization.louvain,
    ):
        G = self.copy()
        G.remove_nodes_from(tuple(nx.isolates(G)))

        if G.number_of_nodes() == 0:
            return []

        communities = algorithm(G, weight=weight)

        while (
            combine_sizes([len(community) for community in communities]) > module_size
        ):
            max_community_size = max(len(community) for community in communities)

            indices = [
                communities.index(community)
                for community in communities
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
                        communities[i : i + 1] = subdivided_community

                        indices[indices.index(i) + 1 :] = [
                            k + len(subdivided_community) - 1
                            for k in indices[indices.index(i) + 1 :]
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
        combine_sites=statistics.mean,
    ):
        proteins = []
        for protein in self:
            sites = [
                self.nodes[protein][change]
                for change in self.nodes[protein]
                if len(change.split(" ")) == 3
                and change.split(" ")[0] == str(time)
                and change.split(" ")[1] == ptm
            ]

            if sites and change_filter(combine_sites(sites)):
                proteins.append(protein)

        return proteins

    def get_module_change_enrichment(
        self,
        p=0.05,
        changes=(-1.0, 1.0),
        get_range=lambda time, ptm, changes, combine_sites: changes,
        combine_sites=statistics.mean,
        module_size=35,
        combine_sizes=statistics.mean,
        weight="weight",
        test="two-sided",
        algorithm=modularization.louvain,
    ):
        modules = {
            i: module
            for i, module in enumerate(
                self.get_modules(
                    module_size, combine_sizes, weight=weight, algorithm=algorithm
                )
            )
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
                        )
                    )

                    k = [
                        len(
                            self.subgraph(modules[i]).get_proteins(
                                time,
                                ptm,
                                lambda change: change > mid_range[1],
                                combine_sites=combine_sites,
                            )
                        )
                        for i in modules
                    ]

                elif test == "one-sided negative":
                    n = len(
                        self.get_proteins(
                            time,
                            ptm,
                            lambda change: change < mid_range[0],
                            combine_sites=combine_sites,
                        )
                    )

                    k = [
                        len(
                            self.subgraph(modules[i]).get_proteins(
                                time,
                                ptm,
                                lambda change: change < mid_range[0],
                                combine_sites=combine_sites,
                            )
                        )
                        for i in modules
                    ]

                elif test == "two-sided":
                    n = len(
                        self.get_proteins(
                            time,
                            ptm,
                            lambda change: change < mid_range[0],
                            combine_sites=combine_sites,
                        )
                    ) + len(
                        self.get_proteins(
                            time,
                            ptm,
                            lambda change: change > mid_range[1],
                            combine_sites=combine_sites,
                        )
                    )

                    k = [
                        len(
                            self.subgraph(modules[i]).get_proteins(
                                time,
                                ptm,
                                lambda change: change < mid_range[0],
                                combine_sites=combine_sites,
                            )
                        )
                        + len(
                            self.subgraph(modules[i]).get_proteins(
                                time,
                                ptm,
                                lambda change: change > mid_range[1],
                                combine_sites=combine_sites,
                            )
                        )
                        for i in modules
                    ]

                else:
                    n = 0
                    k = [0 for _ in modules]

                p_values = {
                    i: scipy.stats.hypergeom.sf(k[i] - 1, M, n, N[i])
                    for i in range(len(modules))
                }

                for module, p_value in correction.benjamini_hochberg(p_values).items():
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
                node for node in self if node.split("-")[0] == protein.split("-")[0]
            }
        else:
            nodes = {protein.split("-")[0]}

        for _ in range(k):
            nodes.update(set.union(*(set(self.neighbors(node)) for node in nodes)))

        return self.subgraph(nodes)
