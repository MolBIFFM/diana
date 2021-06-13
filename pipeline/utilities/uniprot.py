from ..utilities import fetch
from ..configuration import data


def swissprot():
    accessions, entry_gene_name, entry_protein_name = [], {}, {}
    rec_name, taxon_identifier = False, 0
    for line in fetch.txt(data.UNIPROT_SWISSPROT):
        if not line.strip():
            continue

        if line.split(maxsplit=1)[0] == "AC":
            if len(line.split(maxsplit=1)) == 1:
                continue

            accessions.extend(
                line.split(maxsplit=1)[1].rstrip(";").split("; "))

        elif line.split(maxsplit=1)[0] == "GN":
            if len(line.split(maxsplit=1)) == 1:
                continue

            for entry in line.split(maxsplit=1)[1].rstrip(";").split("; "):
                if "=" in entry:
                    entry_gene_name[entry.split("=")[0]] = (
                        entry.split("=")[1].split("{")[0].rstrip())

        elif line.split(maxsplit=1)[0] == "DE":
            if len(line.split(maxsplit=1)) == 1:
                continue

            if line.split(maxsplit=1)[1].split(":", 1)[0] == "RecName":
                entries = (line.split(maxsplit=1)[1].split(
                    ":", 1)[1].lstrip().rstrip(";").split("; "))
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
                            entry.split("=")[1].split("{")[0].rstrip())

        elif line.split(maxsplit=1)[0] == "OX":
            if len(line.split(maxsplit=1)) == 1:
                continue

            if line.split(
                    maxsplit=1)[1].split(";")[0].split("=")[0] == "NCBI_TaxID":
                if (line.split(maxsplit=1)[1].split(";")[0].split("=")
                    [1].split("{")[0].isnumeric()):
                    taxon_identifier = int(
                        line.split(maxsplit=1)[1].split(";")[0].split("=")
                        [1].split("{")[0])

        elif line == "//":
            yield (
                accessions,
                entry_gene_name.get("Name", "NA"),
                entry_protein_name.get("Full", "NA"),
                taxon_identifier,
            )

            accessions.clear()
            entry_gene_name.clear()
            entry_protein_name.clear()
            rec_name, taxon_identifier = False, 0
