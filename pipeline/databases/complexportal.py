import itertools

from uniprot import uniprot
from download import download

COMPLEX_PORTAL = "http://ftp.ebi.ac.uk/pub/databases/intact/complex/current/complextab/{taxon_identifier}.tsv"


def add_proteins(network, evidence_code=[], taxon_identifier=9606):
    primary_accessions = uniprot.get_primary_accession(taxon_identifier)

    nodes_to_add = set()
    for row in download.tabular_txt(
            COMPLEX_PORTAL.format(taxon_identifier=taxon_identifier),
            delimiter="\t",
            header=0,
            usecols=[
                "Taxonomy identifier",
                "Evidence Code",
                "Expanded participant list",
            ],
    ):
        if (row["Taxonomy identifier"] == taxon_identifier
                and (not evidence_code
                     or row["Evidence Code"][row["Evidence Code"].find("(") +
                                             1:row["Evidence Code"].find(")")]
                     in evidence_code)):
            for interactor_a, interactor_b in itertools.combinations(
                (interactor[:interactor.find("(")]
                 for interactor in row["Expanded participant list"].split("|")
                 ), 2):

                if "-" in interactor_a and not interactor_a.split(
                        "-")[1].isnumeric():
                    interactor_a = interactor_a.split("-")[0]

                if "-" in interactor_b and not interactor_b.split(
                        "-")[1].isnumeric():
                    interactor_b = interactor_b.split("-")[0]

                for primary_interactor_a in primary_accessions.get(
                        interactor_a, {interactor_a}):
                    for primary_interactor_b in primary_accessions.get(
                            interactor_b, {interactor_b}):
                        if (primary_interactor_a in network
                                and primary_interactor_b not in network):
                            nodes_to_add.add(primary_interactor_b)

                        elif (primary_interactor_a not in network
                              and primary_interactor_b in network):
                            nodes_to_add.add(primary_interactor_a)

    network.add_nodes_from(nodes_to_add)


def add_interactions(network, evidence_code=[], taxon_identifier=9606):
    primary_accessions = uniprot.get_primary_accession(taxon_identifier,
                                                       network)

    for row in download.tabular_txt(
            COMPLEX_PORTAL.format(taxon_identifier=taxon_identifier),
            delimiter="\t",
            header=0,
            usecols=[
                "Taxonomy identifier",
                "Evidence Code",
                "Expanded participant list",
            ],
    ):
        if (row["Taxonomy identifier"] == taxon_identifier
                and (not evidence_code
                     or row["Evidence Code"][row["Evidence Code"].find("(") +
                                             1:row["Evidence Code"].find(")")]
                     in evidence_code)):
            for interactor_a, interactor_b in itertools.combinations(
                (interactor[:interactor.find("(")]
                 for interactor in row["Expanded participant list"].split("|")
                 ), 2):

                if "-" in interactor_a and not interactor_a.split(
                        "-")[1].isnumeric():
                    interactor_a = interactor_a.split("-")[0]

                if "-" in interactor_b and not interactor_b.split(
                        "-")[1].isnumeric():
                    interactor_b = interactor_b.split("-")[0]

                for primary_interactor_a in primary_accessions.get(
                        interactor_a, {interactor_a}):
                    for primary_interactor_b in primary_accessions.get(
                            interactor_b, {interactor_b}):
                        if (primary_interactor_a in network
                                and primary_interactor_b in network and
                                primary_interactor_a != primary_interactor_b):
                            network.add_edge(
                                primary_interactor_a,
                                primary_interactor_b,
                            )
                            network.edges[
                                primary_interactor_a,
                                primary_interactor_b, ]["ComplexPortal"] = 1.0
