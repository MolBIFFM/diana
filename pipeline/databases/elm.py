from uniprot import uniprot
from download import download

ELM = "http://elm.eu.org/interactions/as_tsv?taxon={organism}"

ORGANISM = {
    9606: "homo%20sapiens",
}


def add_proteins(
    network,
    taxon_identifier=9606,
):
    primary_accession = uniprot.get_primary_accession(taxon_identifier)

    nodes_to_add = set()
    for row in download.tabular_txt(
            ELM.format(organism=ORGANISM[taxon_identifier]),
            delimiter="\t",
            header=0,
            usecols=[
                "interactorElm", "interactorDomain", "taxonomyElm",
                "taxonomyDomain"
            ],
    ):

        if (int(row["taxonomyElm"][:row["taxonomyElm"].find("(")].strip("\""))
                == taxon_identifier
                and int(row["taxonomyDomain"][:row["taxonomyDomain"].find("(")]
                        .strip("\"")) == taxon_identifier):

            interactor_a = row["interactorElm"]
            if "-" in interactor_a and not interactor_a.split(
                    "-")[1].isnumeric():
                interactor_a = interactor_a.split("-")[0]

            for primary_interactor_a in primary_accession.get(
                    interactor_a, {interactor_a}):

                interactor_b = row["interactorDomain"]
                if "-" in interactor_b and not interactor_b.split(
                        "-")[1].isnumeric():
                    interactor_b = interactor_b.split("-")[0]

                for primary_interactor_b in primary_accession.get(
                        interactor_b, {interactor_b}):

                    if (primary_interactor_a in network
                            and primary_interactor_b not in network):
                        nodes_to_add.add(primary_interactor_b)

                    elif (primary_interactor_a not in network
                          and primary_interactor_b in network):
                        nodes_to_add.add(primary_interactor_a)

    network.add_nodes_from(nodes_to_add)


def add_interactions(
    network,
    taxon_identifier=9606,
):
    primary_accession = uniprot.get_primary_accession(taxon_identifier,
                                                      network)

    for row in download.tabular_txt(
            ELM.format(organism=ORGANISM[taxon_identifier]),
            delimiter="\t",
            header=0,
            usecols=[
                "interactorElm", "interactorDomain", "taxonomyElm",
                "taxonomyDomain"
            ],
    ):

        if (int(row["taxonomyElm"][:row["taxonomyElm"].find("(")].strip("\""))
                == taxon_identifier
                and int(row["taxonomyDomain"][:row["taxonomyDomain"].find("(")]
                        .strip("\"")) == taxon_identifier):

            interactor_a = row["interactorElm"]
            if "-" in interactor_a and not interactor_a.split(
                    "-")[1].isnumeric():
                interactor_a = interactor_a.split("-")[0]

            for primary_interactor_a in primary_accession.get(
                    interactor_a, {interactor_a}):

                interactor_b = row["interactorDomain"]
                if "-" in interactor_b and not interactor_b.split(
                        "-")[1].isnumeric():
                    interactor_b = interactor_b.split("-")[0]

                for primary_interactor_b in primary_accession.get(
                        interactor_b, {interactor_b}):

                    if (primary_interactor_a in network
                            and primary_interactor_b in network
                            and primary_interactor_a != primary_interactor_b):
                        network.add_edge(primary_interactor_a,
                                         primary_interactor_b)
                        network.edges[primary_interactor_a,
                                      primary_interactor_b]["ELM"] = 1.0
