"""The interface for the CORUM database."""
import re
from typing import Container, Generator, Optional

from access import iterate

from databases import uniprot

ORGANISM = {"data": {9606: "Homo sapiens"}}


def get_protein_protein_interactions(
    purification_methods: Optional[Container[str]] = None,
    taxonomy_identifier: int = 9606,
) -> Generator[tuple[str, str, float], None, None]:
    """
    Yields protein-protein interactions from CORUM.

    Args:
        purification_methods: The accepted PSI-MI identifiers or terms for the 
            protein complex purification method. If none are specified, any are 
            accepted.
        taxonomy_identifier: The taxonomy identifier.

    Yields:
        Pairs of interacting proteins.
    """
    primary_accession = uniprot.get_primary_accession(taxonomy_identifier)

    for row in iterate.tabular_txt(
            "http://mips.helmholtz-muenchen.de/corum/download/"
            "allComplexes.txt.zip",
            file_from_zip_archive=re.compile(r"allComplexes\.txt\.zip"),
            delimiter="\t",
            header=0,
            usecols=[
                "subunits(UniProt IDs)", "Protein complex purification method",
                "SWISSPROT organism"
            ]):
        if ((not purification_methods or any(
                purification_method.split("-")
            [0].strip() in purification_methods or purification_method.split(
                "-")[1].strip() in purification_methods for purification_method
                in row["Protein complex purification method"].split(";"))) and
                all(organism[:organism.find("(")].rstrip() == ORGANISM["data"]
                    [taxonomy_identifier]
                    for organism in row["SWISSPROT organism"].split(";"))):

            subunits = [
                subunit.split("-")[0] if "-" in subunit and
                not subunit.split("-")[1].isnumeric() else subunit
                for subunit in row["subunits(UniProt IDs)"].split(";")
            ]

            for a in range(len(subunits) - 1):
                for primary_interactor_a in primary_accession.get(
                        subunits[a].split("-")[0], {subunits[a]}):
                    for b in range(a + 1, len(subunits)):
                        for primary_interactor_b in primary_accession.get(
                                subunits[b].split("-")[0], {subunits[b]}):
                            yield (primary_interactor_a, primary_interactor_b)
