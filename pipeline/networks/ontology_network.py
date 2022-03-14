import networkx as nx
from databases import gene_ontology
from enrichment import correction, test


def get_ontology_network(protein_protein_interaction_network,
                         namespaces=("cellular_compartment",
                                     "molecular_function",
                                     "biological_process"),
                         test=test.hypergeometric,
                         correction=correction.benjamini_hochberg,
                         taxon_identifier=9606):
    network = nx.DiGraph()
    for term in gene_ontology.get_ontology(namespaces):
        network.add_node(term["id"])
        network.nodes[term["id"]]["term"] = term["name"]
        network.nodes[term["id"]]["namespace"] = term["namespace"]

        for parent in term.get("is_a", []):
            network.add_edge(term["id"], parent)

    annotation = {term: set() for term in network}
    for protein, term in gene_ontology.get_annotation(taxon_identifier):
        if term in annotation:
            annotation[term].add(protein)

    annotation = {
        term: annotation[term]
        for term in annotation if annotation[term]
    }
    network.remove_nodes_from(term for term in network
                              if term not in annotation)

    annotated_proteins = set.union(*annotation.values())

    intersection = {
        term: len(annotation[term].intersection(
            protein_protein_interaction_network.nodes()))
        for term in annotation
    }

    p_value = correction({
        term: test(
            intersection[term], len(annotated_proteins), len(annotation[term]),
            len(
                annotated_proteins.intersection(
                    protein_protein_interaction_network.nodes())))
        for term in network
    })

    for term in network:
        network.nodes[term]["annotated proteins"] = len(annotation[term])
        network.nodes[term]["network proteins"] = intersection[term]
        network.nodes[term]["p-value"] = p_value[term]

    return network


def get_term_sizes(network):
    return {
        term: network.nodes[term]["annotated proteins"]
        for term in network
    }


def export(network, basename, suffix=""):
    nx.write_graphml_xml(network,
                         "{0}{1}.graphml".format(basename, suffix),
                         named_key_ids=True,
                         infer_numeric_types=True)
