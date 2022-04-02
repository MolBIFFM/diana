"""The interface for UniProt."""
from optparse import Option
from typing import Generator, Optional

from access import iterate


def get_swissprot_entries(
    taxonomy_identifier: Optional[int] = None
) -> Generator[tuple[tuple[str], str, str], None, None]:
    """
    Yields SwissProt entries.

    Args:
        taxonomy_identifier: The taxonomy identifier.

    Yields:
        An entries UniProt accessions, gene name and protein name.
    """
    accessions, entry_gene_name, entry_protein_name = [], {}, {}
    rec_name, tax_id = False, 0
    for line in iterate.txt(
            "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.dat.gz"
    ):
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
                if (line.split(maxsplit=1)[1].split(";")[0].split("=")[1].split(
                        "{")[0].isnumeric()):
                    tax_id = int(
                        line.split(maxsplit=1)[1].split(";")[0].split("=")
                        [1].split("{")[0])

        elif line == "//":
            if not taxonomy_identifier or tax_id == taxonomy_identifier:
                yield (tuple(accessions), entry_gene_name.get("Name", "NA"),
                       entry_protein_name.get("Full", "NA"))

            accessions.clear()
            entry_gene_name.clear()
            entry_protein_name.clear()
            rec_name, tax_id = False, 0


def get_primary_accession(
        taxonomy_identifier: Optional[int] = None) -> dict[str, set[str]]:
    """
    Returns a map of primary UniProt accessions.

    Args:
        taxonomy_identifier: The taxonomy identifier.

    Returns:
        A map of any secondary UniProt accessions to its primary equivalents.
    """
    primary_accession = {}
    for accessions, _, _ in get_swissprot_entries(taxonomy_identifier):
        for i, accession in enumerate(accessions):
            if i > 0:
                if accession not in primary_accession:
                    primary_accession[accession] = set()
                primary_accession[accession].add(accessions[0])

    return primary_accession
