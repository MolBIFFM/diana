"""The interface for the BioGRID database."""
import re
from typing import Container, Generator, Optional

from access import iterate

from databases import uniprot

ORGANISM = {"files": {9606: "Homo_sapiens"}}


def get_protein_protein_interactions(
    experimental_system: Optional[Container[str]] = None,
    experimental_system_type: Optional[Container[str]] = None,
    interaction_throughput: Optional[Container[str]] = None,
    multi_validated_physical: bool = False,
    taxonomy_identifier: int = 9606,
    version: Optional[tuple[int, int, int]] = None
) -> Generator[tuple[str, str], None, None]:
    """
    Yields protein-protein interactions from BioGRID.

    Args:
        experimental_system: The accepted experimental evidence codes. If none
            are specified, any is accepted.
        experimental_system_type: The accepted categories of experimental
            evidence. If none are specified, any is accepted.
        interaction_throughput:  The accepted levels of interaction throughput.
            If none are specified, any is accepted.
        multi-validated physical: If True, yield only multi-validated physical
            interactions.
        taxonomy_identifier: The taxonomy identifier.
        version: The version of the BioGRID database, if not passed, the latest.

    Yields:
        Pairs of interacting proteins.
    """
    primary_accession = uniprot.get_primary_accession(taxonomy_identifier)

    for row in iterate.tabular_txt(
        ("https://downloads.thebiogrid.org/Download/BioGRID/Release-Archive/"
         f"BIOGRID-{version[0]}.{version[1]}.{version[2]}/"
         f"BIOGRID-MV-Physical-{version[0]}.{version[1]}.{version[2]}.tab3.zip"
         if multi_validated_physical else
         "https://downloads.thebiogrid.org/Download/BioGRID/Release-Archive/"
         f"BIOGRID-{version[0]}.{version[1]}.{version[2]}/"
         f"BIOGRID-ORGANISM-{version[0]}.{version[1]}.{version[2]}.tab3.zip")
            if version else
        ("https://downloads.thebiogrid.org/Download/BioGRID/Latest-Release/"
         "BIOGRID-MV-Physical-LATEST.tab3.zip" if multi_validated_physical else
         "https://downloads.thebiogrid.org/Download/BioGRID/Latest-Release/"
         "BIOGRID-ORGANISM-LATEST.tab3.zip"),
            file_from_zip_archive=re.compile(
                r"BIOGRID-MV-Physical-[0-9]\.[0-9]\.[0-9]{3}\.tab3\.txt"
                if multi_validated_physical else
                f"BIOGRID-ORGANISM-{ORGANISM['files'][taxonomy_identifier]}-"
                r"[0-9]\.[0-9]\.[0-9][0-9][0-9]\.tab3\.txt"),
            delimiter="\t",
            header=0,
            usecols=[
                "Experimental System", "Experimental System Type",
                "Organism ID Interactor A", "Organism ID Interactor B",
                "Throughput", "SWISS-PROT Accessions Interactor A",
                "SWISS-PROT Accessions Interactor B"
            ],
    ):
        if (row["SWISS-PROT Accessions Interactor A"] == "-" or
                row["SWISS-PROT Accessions Interactor B"] == "-"):
            continue

        if ((not experimental_system or
             row["Experimental System"] in experimental_system) and
            (not experimental_system_type or
             row["Experimental System Type"] in experimental_system_type) and
            (row["Organism ID Interactor A"] == taxonomy_identifier and
             row["Organism ID Interactor B"] == taxonomy_identifier) and
            (not interaction_throughput or
             any(it in interaction_throughput
                 for it in row["Throughput"].split("|")))):
            for interactor_a in row["SWISS-PROT Accessions Interactor A"].split(
                    "|"):
                if "-" in interactor_a and not interactor_a.split(
                        "-")[1].isnumeric():
                    interactor_a = interactor_a.split("-")[0]
                for interactor_b in row[
                        "SWISS-PROT Accessions Interactor B"].split("|"):
                    if "-" in interactor_b and not interactor_b.split(
                            "-")[1].isnumeric():
                        interactor_b = interactor_b.split("-")[0]
                    for primary_interactor_a in primary_accession.get(
                            interactor_a.split("-")[0], {interactor_a}):
                        for primary_interactor_b in primary_accession.get(
                                interactor_b.split("-")[0], {interactor_b}):
                            yield (primary_interactor_a, primary_interactor_b)
