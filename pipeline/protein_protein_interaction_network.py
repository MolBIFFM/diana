import bisect
import itertools
import json
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

    def add_genes_from_spreadsheet(
        self,
        file_name,
        gene_accession_column,
        gene_accession_format,
        sheet_name=0,
        header=0,
        organism=9606,
    ):
        genes = []
        for _, row in pd.read_excel(
            file_name,
            sheet_name=sheet_name,
            header=header,
            usecols=[gene_accession_column],
            dtype={gene_accession_column: str},
        ).iterrows():

            if pd.isna(row[gene_accession_column]):
                continue

            genes.extend(
                [
                    str(gene_accession)
                    for gene_accession in gene_accession_format(
                        row[gene_accession_column]
                    )
                ]
            )

        proteins, gene_name, protein_name = [], {}, {}

        accessions, entry_gene_name, entry_protein_name = [], {}, {}
        rec_name, ncbi_tax_id = False, 0
        for line in fetch.txt(data.UNIPROT_SWISS_PROT):
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
                        ncbi_tax_id = int(
                            line.split(maxsplit=1)[1]
                            .split(";")[0]
                            .split("=")[1]
                            .split("{")[0]
                        )

            elif line == "//":
                if ncbi_tax_id == organism and entry_gene_name.get("Name") in genes:
                    proteins.append(accessions[0])
                    gene_name[accessions[0]] = entry_gene_name["Name"]
                    protein_name[accessions[0]] = entry_protein_name.get("Full", "NA")

                accessions.clear()
                entry_gene_name.clear()
                entry_protein_name.clear()

                rec_name, ncbi_tax_id = False, 0

        for protein in proteins:
            self.add_node(protein)
            self.nodes[protein]["gene name"] = gene_name[protein]
            self.nodes[protein]["protein name"] = protein_name[protein]

            yield (
                self.nodes[protein]["gene name"],
                protein,
                self.nodes[protein]["protein name"],
            )

    def add_proteins_from_spreadsheet(
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
        merge_replicates=statistics.mean,
        convert_measurement=math.log2,
        organism=9606,
    ):
        proteins = {}
        for _, row in pd.read_excel(
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
        ).iterrows():

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
                            protein, isoform = protein_accession, None

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
                        protein, isoform = protein_accession, None

                    if protein not in proteins:
                        proteins[protein] = {}

                    if isoform not in proteins[protein]:
                        proteins[protein][isoform] = []

        swissprot_proteins, primary_accession, gene_name, protein_name = {}, {}, {}, {}

        accessions, entry_gene_name, entry_protein_name = [], {}, {}
        rec_name, ncbi_tax_id = False, 0
        for line in fetch.txt(data.UNIPROT_SWISS_PROT):
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
                        ncbi_tax_id = int(
                            line.split(maxsplit=1)[1]
                            .split(";")[0]
                            .split("=")[1]
                            .split("{")[0]
                        )

            elif line == "//":
                if ncbi_tax_id == organism:
                    for i, accession in enumerate(accessions):
                        if accession in proteins:
                            swissprot_proteins[accession] = proteins[accession]
                            gene_name[accession] = entry_gene_name.get("Name", "NA")
                            protein_name[accession] = entry_protein_name.get(
                                "Full", "NA"
                            )
                            if i > 0:
                                primary_accession[accession] = accessions[0]

                accessions.clear()
                entry_gene_name.clear()
                entry_protein_name.clear()

                rec_name, ncbi_tax_id = False, 0

        proteins = swissprot_proteins

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
                                primary_position
                                if primary_position
                                else secondary_position
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
                if isoform:
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

                proteins[protein][isoform] = sorted(
                    sorted(
                        proteins[protein][isoform],
                        key=lambda item: item[1],
                        reverse=True,
                    )[:num_sites]
                )

                for i in range(len(proteins[protein][isoform])):
                    if isoform:
                        self.nodes["-".join([protein, isoform])][
                            "{} {} {}".format(time, ptm, i + 1)
                        ] = proteins[protein][isoform][i][1]
                    else:
                        self.nodes[protein][
                            "{} {} {}".format(time, ptm, i + 1)
                        ] = proteins[protein][isoform][i][1]

                yield (
                    self.nodes["-".join([protein, isoform]) if isoform else protein][
                        "gene name"
                    ],
                    "-".join([protein, isoform]) if isoform else protein,
                    self.nodes["-".join([protein, isoform]) if isoform else protein][
                        "protein name"
                    ],
                )

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

    def get_changes(self, time, ptm, merge_sites=statistics.mean):
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
                changes.append(merge_sites(sites))

        return changes

    def get_z_score_range(
        self,
        time,
        post_translational_modification,
        thresholds,
        merge_sites=statistics.mean,
    ):
        changes = self.get_changes(time, post_translational_modification, merge_sites)
        mean = statistics.mean(changes)
        stdev = statistics.stdev(changes, xbar=mean)
        return (thresholds[0] * stdev + mean, thresholds[1] * stdev + mean)

    def get_z_score(
        self,
        time,
        post_translational_modification,
        change,
        merge_sites=statistics.mean,
    ):
        changes = self.get_changes(time, post_translational_modification, merge_sites)
        mean = statistics.mean(changes)
        stdev = statistics.stdev(changes, xbar=mean)
        return (change - mean) / stdev

    def get_proportion_range(
        self,
        time,
        post_translational_modification,
        thresholds,
        merge_sites=statistics.mean,
    ):
        changes = sorted(
            self.get_changes(time, post_translational_modification, merge_sites)
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
        merge_sites=statistics.mean,
    ):
        changes = self.get_changes(time, post_translational_modification, merge_sites)
        return (
            min(
                len([c for c in changes if c <= change]),
                len([c for c in changes if c >= change]),
            )
            / len(changes)
        )

    def set_changes(
        self,
        merge_sites=statistics.mean,
        changes=(-1.0, 1.0),
        get_range=lambda time, post_translational_modification, changes, merge_sites: changes,
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
                    merge_sites,
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
                else:
                    self.nodes[protein]["change {}".format(time)] = "mid"

    def get_databases(self):
        return sorted(
            set(
                database for edge in self.edges for database in self.edges[edge]
            ).intersection({"BioGRID", "IntAct", "STRING", "Reactome", "CORUM"})
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

    def add_interactions_from_biogrid(
        self, experimental_system=[], organism=9606, multi_validated_physical=False
    ):
        uniprot = {}
        for row in fetch.tabular_txt(
            data.UNIPROT_ID_MAP,
            delimiter="\t",
            usecols=[0, 1, 2],
        ):
            if row[1] == "BioGRID" and row[0] in self.nodes:
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
                "ORGANISM_OFFICIAL_NAME",
            ],
        ):
            if (
                row["IDENTIFIER_TYPE"] in ("UNIPROT-ACCESSION", "UNIPROT-ISOFORM")
                and row["ORGANISM_OFFICIAL_NAME"] == "Homo sapiens"
                and row["IDENTIFIER_VALUE"] in self.nodes
            ):
                uniprot[int(row["BIOGRID_ID"])] = row["IDENTIFIER_VALUE"]

        for row in fetch.tabular_txt(
            data.BIOGRID_MV_PHYSICAL_ZIP_ARCHIVE
            if multi_validated_physical
            else data.BIOGRID_ZIP_ARCHIVE,
            file=data.BIOGRID_MV_PHYSICAL if multi_validated_physical else data.BIOGRID,
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
                and row["Organism ID Interactor A"] == organism
                and row["Organism ID Interactor B"] == organism
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

    def add_interactions_from_corum(
        self, protein_complex_purification_method=[], organism="Human"
    ):
        for row in fetch.tabular_txt(
            data.CORUM_ZIP_ARCHIVE,
            file=data.CORUM,
            delimiter="\t",
            header=0,
            usecols=[
                "Organism",
                "subunits(UniProt IDs)",
                "Protein complex purification method",
            ],
        ):
            if row["Organism"] == organism and (
                not protein_complex_purification_method
                or any(
                    method in protein_complex_purification_method
                    for method in [
                        entry.split("-")[1].lstrip()
                        for entry in row["Protein complex purification method"].split(
                            ";"
                        )
                    ]
                )
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

                        yield (
                            interactor_a,
                            interactor_b,
                            self.edges[
                                interactor_a,
                                interactor_b,
                            ]["CORUM"],
                        )

    def add_interactions_from_intact(
        self,
        interaction_detection_methods=[],
        interaction_types=[],
        mi_score=0.0,
        organism=9606,
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

            if not (
                interactor_a := mitab.get_identifier_from_namespace(
                    row["#ID(s) interactor A"], "uniprotkb"
                )
            ):
                if not (
                    interactor_a := mitab.get_identifier_from_namespace(
                        row["Alt. ID(s) interactor A"], "uniprotkb"
                    )
                ):
                    continue
            if interactor_a not in self:
                continue

            if not (
                interactor_b := mitab.get_identifier_from_namespace(
                    row["ID(s) interactor B"], "uniprotkb"
                )
            ):
                if not (
                    interactor_b := mitab.get_identifier_from_namespace(
                        row["Alt. ID(s) interactor B"], "uniprotkb"
                    )
                ):
                    continue
            if interactor_b not in self:
                continue

            if interactor_a == interactor_b:
                continue

            if not (
                mitab.namespace_has_identifier(
                    row["Taxid interactor A"], "taxid", organism
                )
                and mitab.namespace_has_identifier(
                    row["Taxid interactor B"], "taxid", organism
                )
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

            if score := mitab.get_identifier_from_namespace(
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

    def add_interactions_from_reactome(
        self,
        interaction_detection_methods=[],
        interaction_types=[],
        reactome_score=0.0,
        organism=9606,
    ):
        for row in fetch.tabular_txt(
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
                interactor_a := mitab.get_identifier_from_namespace(
                    row["#ID(s) interactor A"], "uniprotkb"
                )
            ):
                if not (
                    interactor_a := mitab.get_identifier_from_namespace(
                        row["Alt. ID(s) interactor A"], "uniprotkb"
                    )
                ):
                    continue
            if interactor_a not in self:
                continue

            if not (
                interactor_b := mitab.get_identifier_from_namespace(
                    row["ID(s) interactor B"], "uniprotkb"
                )
            ):
                if not (
                    interactor_b := mitab.get_identifier_from_namespace(
                        row["Alt. ID(s) interactor B"], "uniprotkb"
                    )
                ):
                    continue
            if interactor_b not in self:
                continue

            if interactor_a == interactor_b:
                continue

            if not (
                mitab.namespace_has_identifier(
                    row["Taxid interactor A"], "taxid", organism
                )
                and mitab.namespace_has_identifier(
                    row["Taxid interactor B"], "taxid", organism
                )
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

            if score := mitab.get_identifier_from_namespace(
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
        organism=9606,
        physical=False,
    ):

        uniprot = {}
        for row in fetch.tabular_txt(
            data.UNIPROT_ID_MAP,
            delimiter="\t",
            usecols=[0, 1, 2],
        ):
            if row[1] == "STRING" and row[0] in self.nodes:
                uniprot[row[2]] = row[0]

        for row in fetch.tabular_txt(data.STRING_ID_MAP, usecols=[1, 2]):
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

        for row in fetch.tabular_txt(
            data.STRING_PHYSICAL.format(organism=organism)
            if physical
            else data.STRING.format(organism=organism),
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

    def get_modules(
        self,
        module_size=35,
        merge_sizes=statistics.mean,
        weight="weight",
        algorithm=modularization.louvain,
    ):
        self.remove_nodes_from(list(nx.isolates(self)))

        if self.number_of_nodes() == 0:
            return []

        communities = [
            community
            for community in list(nx.algorithms.components.connected_components(self))
        ]

        while merge_sizes([len(community) for community in communities]) > module_size:
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
                        (self.subgraph(communities[j]) for j in indices),
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
        change_filter=lambda merged_sites: bool(merged_sites),
        merge_sites=statistics.mean,
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

            if sites and change_filter(merge_sites(sites)):
                proteins.append(protein)

        return proteins

    def get_module_change_enrichment(
        self,
        p=0.05,
        changes=(-1.0, 1.0),
        get_range=lambda time, ptm, changes, merge_sites: changes,
        merge_sites=statistics.mean,
        module_size=35,
        merge_sizes=statistics.mean,
        weight="weight",
        test="two-sided",
        algorithm=modularization.louvain,
    ):
        modules = {
            i: module
            for i, module in enumerate(
                self.get_modules(
                    module_size, merge_sizes, weight=weight, algorithm=algorithm
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
                    merge_sites,
                )

                if test == "one-sided positive":
                    n = len(
                        self.get_proteins(
                            time,
                            ptm,
                            lambda change: change > mid_range[1],
                            merge_sites=merge_sites,
                        )
                    )

                    k = [
                        len(
                            self.subgraph(modules[i]).get_proteins(
                                time,
                                ptm,
                                lambda change: change > mid_range[1],
                                merge_sites=merge_sites,
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
                            merge_sites=merge_sites,
                        )
                    )

                    k = [
                        len(
                            self.subgraph(modules[i]).get_proteins(
                                time,
                                ptm,
                                lambda change: change < mid_range[0],
                                merge_sites=merge_sites,
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
                            merge_sites=merge_sites,
                        )
                    ) + len(
                        self.get_proteins(
                            time,
                            ptm,
                            lambda change: change > mid_range[1],
                            merge_sites=merge_sites,
                        )
                    )

                    k = [
                        len(
                            self.subgraph(modules[i]).get_proteins(
                                time,
                                ptm,
                                lambda change: change < mid_range[0],
                                merge_sites=merge_sites,
                            )
                        )
                        + len(
                            self.subgraph(modules[i]).get_proteins(
                                time,
                                ptm,
                                lambda change: change > mid_range[1],
                                merge_sites=merge_sites,
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
