from databases import uniprot
from fetch import fetch
from enrichment import test, correction

ONTOLOGY = "http://purl.obolibrary.org/obo/go.obo"
ANNOTATION = "http://geneontology.org/gene-associations/goa_{organism}.gaf.gz"
ANNOTATION_ISOFORM = "http://geneontology.org/gene-associations/goa_{organism}_isoform.gaf.gz"

ORGANISM = {"file": {9606: "human"}}


def get_ontology(namespaces=("biological_process", "cellular_compartment",
                             "molecular_function")):
    term = {}
    for line in fetch.txt(ONTOLOGY):
        if any(
                line.startswith("{}:".format(tag))
                for tag in ("format-version", "data-version", "subsetdef",
                            "synonymtypedef", "default-namespace", "ontology",
                            "property_value")):
            continue
        elif line == "[Term]" or line == "[Typedef]":
            if term.get("namespace") in namespaces:
                yield term
            term = {}
        elif any(
                line.startswith("{}:".format(tag))
                for tag in ("id", "name", "namespace")):
            term[line.split(":", maxsplit=1)[0]] = line.split(
                ":", maxsplit=1)[1].strip()
        elif line.startswith("is_a:"):
            if "is_a" not in term:
                term["is_a"] = []
            term["is_a"].append(
                line.split(":", maxsplit=1)[1].split("!")[0].strip())


def get_annotation(taxon_identifier=9606):
    primary_accession = uniprot.get_primary_accession(taxon_identifier)

    for row in fetch.tabular_txt(
            ANNOTATION.format(organism=ORGANISM["file"][taxon_identifier]),
            skiprows=41,
            delimiter="\t",
            usecols=[0, 1, 4, 12]):
        if row[0] == "UniProtKB" and row[12] == "taxon:{}".format(
                taxon_identifier):
            for protein in primary_accession.get(row[1], {row[1]}):
                yield (protein, row[4])

    for row in fetch.tabular_txt(ANNOTATION_ISOFORM.format(
            organism=ORGANISM["file"][taxon_identifier]),
                                 skiprows=41,
                                 delimiter="\t",
                                 usecols=[0, 4, 12, 16]):
        if row[0] == "UniProtKB" and row[12] == "taxon:{}".format(
                taxon_identifier) and row[16].startswith("UniProtKB:"):
            yield (row[16].split(":")[1], row[4])


def get_enrichment(networks,
                   test=test.hypergeometric,
                   correction=correction.benjamini_hochberg,
                   taxon_identifier=9606,
                   namespaces=("biological_process", "cellular_compartment",
                               "molecular_function")):
    annotation = {}
    for protein, term in get_annotation(taxon_identifier):
        if term not in annotation:
            annotation[term] = set()
        annotation[term].add(protein)

    name = {}
    for term in get_ontology(namespaces):
        if term["id"] in annotation:
            name[term["id"]] = term["name"]

    annotation = {
        term: proteins
        for term, proteins in annotation.items() if proteins and term in name
    }

    annotated_proteins = set.union(*annotation.values())

    intersection = {
        network: {
            term: len(annotation[term].intersection(network.nodes()))
            for term in annotation
        }
        for network in networks
    }

    p_values = correction({
        (network, term):
        test(intersection[network][term], len(annotated_proteins),
             len(annotation[term]),
             len(annotated_proteins.intersection(network.nodes())))
        for term in annotation for network in networks
    })

    return {
        network: {(term, name[term]): p_values[(network, term)]
                  for term in annotation}
        for network in networks
    }
