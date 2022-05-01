"""The interface for the UniProt database."""
from typing import Iterator, Optional

from access import iterate


def get_swiss_prot_entries(
        organism: Optional[int] = None
) -> Iterator[tuple[tuple[str], str, str]]:
    """
    Yields Swiss-Prot entries.

    Args:
        organism: The NCBI taxonomy identifier for the organism of interest. 

    Yields:
        An entries UniProt accessions, gene name and protein name.
    """
    accessions, entry_gene_name, entry_protein_name = [], {}, {}
    rec_name, tax_id = False, 0
    for line in iterate.txt(
            "https://ftp.uniprot.org/pub/databases/uniprot/current_release/"
            "knowledgebase/complete/uniprot_sprot.dat.gz"):
        if not line.strip():
            continue

        if line.split()[0] == "AC":
            if len(line.split()) == 1:
                continue

            accessions.extend(
                line.split(maxsplit=1)[1].rstrip(";").split("; "))

        elif line.split()[0] == "GN":
            if len(line.split()) == 1:
                continue

            for entry in line.split(maxsplit=1)[1].rstrip(";").split("; "):
                if "=" in entry:
                    entry_gene_name[entry.split("=")[0]] = (
                        entry.split("=")[1].split("{")[0].rstrip())

        elif line.split()[0] == "DE":
            if len(line.split()) == 1:
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

        elif line.split()[0] == "OX":
            if len(line.split()) == 1:
                continue

            if line.split(
                    maxsplit=1)[1].split(";")[0].split("=")[0] == "NCBI_TaxID":
                if (line.split(maxsplit=1)[1].split(";")[0].split("=")
                    [1].split("{")[0].isnumeric()):
                    tax_id = int(
                        line.split(maxsplit=1)[1].split(";")[0].split("=")
                        [1].split("{")[0])

        elif line == "//":
            if not organism or tax_id == organism:
                yield (tuple(accessions), entry_gene_name.get("Name", "NA"),
                       entry_protein_name.get("Full", "NA"))

            accessions.clear()
            entry_gene_name.clear()
            entry_protein_name.clear()
            rec_name, tax_id = False, 0


def get_primary_accession(
        organism: Optional[int] = None) -> dict[str, frozenset[str]]:
    """
    Returns a map of primary UniProt accessions.

    Args:
        organism: The NCBI taxonomy identifier for the organism of interest. 

    Returns:
        A map of any secondary UniProt accessions to its primary equivalents.
    """
    primary_accession = {}
    for accessions, _, _ in get_swiss_prot_entries(organism):
        for i, accession in enumerate(accessions):
            if i > 0:
                if accession not in primary_accession:
                    primary_accession[accession] = set()
                primary_accession[accession].add(accessions[0])

    return {
        accession: frozenset(accessions)
        for accession, accessions in primary_accession.items()
    }
