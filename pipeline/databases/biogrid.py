from typing import Generator
from databases import uniprot
from download import download

ORGANISM = {"file": {9606: "Homo_sapiens"}}


def get_protein_protein_interactions(
    experimental_system: list[str] = [],
    experimental_system_type: list[str] = [],
    interaction_throughput: list[str] = [],
    multi_validated_physical: bool = False,
    taxonomy_identifier: int = 9606,
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

    Yields:
        Pairs of interacting proteins.
    """
    primary_accession = uniprot.get_primary_accession(taxonomy_identifier)

    for row in download.tabular_txt(
            "https://downloads.thebiogrid.org/Download/BioGRID/Latest-Release/BIOGRID-MV-Physical-LATEST.tab3.zip"
            if multi_validated_physical else
            "https://downloads.thebiogrid.org/Download/BioGRID/Latest-Release/BIOGRID-ORGANISM-LATEST.tab3.zip",
            file_from_zip_archive=
            r"BIOGRID-MV-Physical-[0-9]\.[0-9]\.[0-9]{3}\.tab3\.txt"
            if multi_validated_physical else
            r"BIOGRID-ORGANISM-{organism}-[0-9]\.[0-9]\.[0-9][0-9][0-9]\.tab3\.txt"
            .format(organism=ORGANISM["file"][taxonomy_identifier]),
            delimiter="\t",
            header=0,
            usecols=[
                "Experimental System", "Experimental System Type", "Throughput",
                "SWISS-PROT Accessions Interactor A",
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
                            interactor_a, {interactor_a}):
                        for primary_interactor_b in primary_accession.get(
                                interactor_b, {interactor_b}):
                            yield (primary_interactor_a, primary_interactor_b)
