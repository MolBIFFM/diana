"""The interface for the CORUM database."""
import re
from typing import (Callable, Container, Hashable, Iterable, Iterator, Mapping,
                    Optional)

import scipy.stats
from access import iterate
from algorithms import correction

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


def get_protein_complexes(
    purification_methods: Optional[Container[str]] = None,
    organism: int = 9606,
) -> Iterator[tuple[int, str, frozenset[str]]]:
    """
    Yields protein complexes from CORUM.

    Args:
        purification_methods: The accepted PSI-MI identifiers or terms for the
            protein complex purification method. If none are specified, any are
            accepted.
        organism: The NCBI taxonomy identifier for the organism of interest.

    Yields:
        Identifier, name and subunits of protein complexes.
    """
    primary_accession = uniprot.get_primary_accession(organism)

    for row in iterate.tabular_txt(
            "https://mips.helmholtz-muenchen.de/corum/download/releases/"
            "current/allComplexes.txt.zip",
            file_from_zip_archive=re.compile(r"allComplexes\.txt\.zip"),
            delimiter="\t",
            header=0,
            usecols=[
                "ComplexID", "ComplexName", "subunits(UniProt IDs)",
                "Protein complex purification method", "SWISSPROT organism"
            ]):
        if ((not purification_methods or any(
                purification_method.split("-")
            [0].strip() in purification_methods or purification_method.split(
                "-")[1].strip() in purification_methods for purification_method
                in row["Protein complex purification method"].split(";"))) and
                all(spo[:spo.find("(")].strip() == ORGANISM["data"][organism]
                    for spo in row["SWISSPROT organism"].split(";"))):

            subunits = [
                subunit.split("-")[0] if "-" in subunit and
                not subunit.split("-")[1].isnumeric() else subunit
                for subunit in row["subunits(UniProt IDs)"].split(";")
            ]

            yield (row["ComplexID"], row["ComplexName"],
                   frozenset.union(
                       *(frozenset(primary_accession.get(subunit, {subunit}))
                         for subunit in subunits)))


def get_enrichment(
    proteins: Iterable[frozenset[str]],
    enrichment_test: Callable[
        [int, int, int, int],
        float] = lambda k, M, n, N: scipy.stats.hypergeom.sf(k - 1, M, n, N),
    multiple_testing_correction: Callable[[dict[Hashable, float]], Mapping[
        Hashable, float]] = correction.benjamini_hochberg,
    purification_methods: Optional[Container[str]] = None,
    organism: int = 9606,
    annotation_as_reference: bool = False
) -> dict[frozenset[str], dict[tuple[int, str], float]]:
    """
    Test sets of proteins for enrichment of CORUM protein complexes.

    Args:
        proteins: The sets of proteins to test.
        enrichment_test: The statistical test used to assess enrichment of
            CORUM protein complexes.
        multiple_testing_correction: The procedure to correct for multiple
            testing of multiple pathways and sets of proteins.
        purification_methods: The accepted PSI-MI identifiers or terms for the
            protein complex purification method. If none are specified, any are
            accepted.
        organism: The NCBI taxonomy identifier for the organism of interest.
        annotation_as_reference: If True, compute enrichment with respect to the
            species-specific set of protein complexes, otherwise with respect to
            the union of the sets of proteins.

    Returns:
        Corrected p-value for the enrichment of each CORUM protein complex by
        each set of proteins.
    """
    annotation, name = {}, {}
    for protein_complex, complex_name, subunits in get_protein_complexes(
            purification_methods, organism):
        if annotation_as_reference or any(
                subunits.intersection(prt) for prt in proteins):
            annotation[protein_complex] = subunits
            name[protein_complex] = complex_name

    annotated_proteins = frozenset.union(*annotation.values())

    annotated_network_proteins = {
        prt: len(annotated_proteins.intersection(prt)) for prt in proteins
    }

    network_intersection = {
        prt: {
            protein_complex: len(subunits.intersection(prt))
            for protein_complex, subunits in annotation.items()
        } for prt in proteins
    }

    p_value = multiple_testing_correction({
        (prt, protein_complex):
        enrichment_test(network_intersection[prt][protein_complex],
                        len(annotated_proteins), len(subunits),
                        annotated_network_proteins[prt])
        for protein_complex, subunits in annotation.items() for prt in proteins
    })

    return {
        prt: {(protein_complex, name[protein_complex]):
              p_value[(prt, protein_complex)]
              for protein_complex in annotation} for prt in proteins
    }
