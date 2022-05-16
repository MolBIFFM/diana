"""The interface for the UniProt database."""
from typing import Iterator

from access import iterate

ORGANISM = {"files": {9606: "human"}}


def get_swiss_prot_entries(
        organism: int = 9606) -> Iterator[tuple[tuple[str, ...], str, str]]:
    """
    Yields Swiss-Prot entries.

    Args:
        organism: The NCBI taxonomy identifier for the organism of interest.

    Yields:
        UniProt accessions, gene name and protein name of each entry.
    """
    accessions = []
    entry_gene_name, entry_protein_name, entry_taxonomy = {}, {}, {}
    rec_name = False

    for line in iterate.txt(
            "https://ftp.uniprot.org/pub/databases/uniprot/current_release/"
            "knowledgebase/taxonomic_divisions/"
            f"uniprot_sprot_{ORGANISM['files'][organism]}.dat.gz"):

        if line.startswith("AC") and len(line.split(maxsplit=1)) == 2:
            accessions.extend(
                line.split(maxsplit=1)[1].rstrip(";").split("; "))

        elif line.startswith("DE") and len(line.split(maxsplit=1)) == 2:
            if line.split(maxsplit=1)[1].split(":", 1)[0] == "RecName":
                entries = line.split(maxsplit=1)[1].split(
                    ":", 1)[1].lstrip().rstrip(";").split("; ")
                rec_name = True

            elif rec_name:
                entries = line.split(maxsplit=1)[1].rstrip(";").split("; ")

            elif any(
                    line.split(maxsplit=1)[1].split(":", 1)[0] == tag
                    for tag in ("AltName", "Contains", "Flags")):
                entries = []
                rec_name = False

            for entry in entries:
                if "=" in entry:
                    entry_protein_name[entry.split("=")[0]] = (
                        entry.split("=")[1].split("{")[0].rstrip())

        elif line.startswith("GN") and len(line.split(maxsplit=1)) == 2:
            for entry in line.split(maxsplit=1)[1].rstrip(";").split("; "):
                if "=" in entry:
                    entry_gene_name[entry.split("=")[0]] = (
                        entry.split("=")[1].split("{")[0].rstrip())

        elif line.startswith("OX") and len(line.split(maxsplit=1)) == 2:
            for entry in line.split(maxsplit=1)[1].rstrip(";").split("; "):
                if "=" in entry:
                    entry_taxonomy[entry.split("=")[0]] = (
                        entry.split("=")[1].split("{")[0].rstrip())

        elif line.startswith("//"):
            if int(entry_taxonomy.get("NCBI_TaxID", 0)) == organism:
                yield (tuple(accessions), entry_gene_name.get("Name", "NA"),
                       entry_protein_name.get("Full", "NA"))

            accessions.clear()
            entry_gene_name.clear()
            entry_protein_name.clear()
            entry_taxonomy.clear()
            rec_name = False


def get_primary_accession(organism: int = 9606) -> dict[str, frozenset[str]]:
    """
    Returns a map of primary UniProt accessions.

    Args:
        organism: The NCBI taxonomy identifier for the organism of interest.

    Returns:
        A map of any secondary UniProt accessions to its primary equivalents.
    """
    primary_accession: dict[str, set[str]] = {}
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
