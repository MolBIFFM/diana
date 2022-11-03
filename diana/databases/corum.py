"""The interface for the CORUM database."""
import re
from typing import Container, Iterator, Optional

from access import iterate

from databases import uniprot

ORGANISM: dict[str, dict[int, str]] = {"data": {9606: "Homo sapiens"}}


def get_protein_interactions(
    purification_methods: Optional[Container[str]] = None,
    organism: int = 9606,
) -> Iterator[tuple[str, str]]:
    """
    Yields protein-protein interactions from CORUM.

    Args:
        purification_methods: The accepted PSI-MI identifiers or terms for the
            protein complex purification method. If none are specified, any are
            accepted.
        organism: The NCBI taxonomy identifier for the organism of interest.

    Yields:
        Pairs of interacting proteins.
    """
    primary_accession = uniprot.get_primary_accession(organism)

    for row in iterate.tabular_txt(
            "https://mips.helmholtz-muenchen.de/corum/download/releases/"
            "current/allComplexes.txt.zip",
            file_from_zip_archive=re.compile(r"^allComplexes\.txt\.zip$"),
            delimiter="\t",
            header=0,
            usecols=[
                "subunits(UniProt IDs)", "Protein complex purification method",
                "SWISSPROT organism"
            ]):
        if (not purification_methods or any(
                purification_method.split("-")
            [0].strip() in purification_methods or purification_method.split(
                "-")[1].strip() in purification_methods for purification_method
                in row["Protein complex purification method"].split(";"))):

            subunits = [
                subunit.split("-")[0] if "-" in subunit and
                not subunit.split("-")[1].isnumeric() else subunit
                for subunit in row["subunits(UniProt IDs)"].split(";")
            ]

            organisms = row["SWISSPROT organism"].split(";")

            for a in range(len(subunits) - 1):
                if organisms[a][:organisms[a].find("(")].strip(
                ) != ORGANISM["data"][organism]:
                    continue
                for primary_interactor_a in primary_accession.get(
                        subunits[a].split("-")[0], {subunits[a]}):
                    for b in range(a + 1, len(subunits)):
                        if organisms[b][:organisms[b].find("(")].strip(
                        ) != ORGANISM["data"][organism]:
                            continue
                        for primary_interactor_b in primary_accession.get(
                                subunits[b].split("-")[0], {subunits[b]}):
                            yield (primary_interactor_a, primary_interactor_b)
