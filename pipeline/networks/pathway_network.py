import networkx as nx
from databases import reactome
from enrichment import correction, test


def get_pathway_network(protein_protein_interaction_network,
                        test=test.hypergeometric,
                        correction=correction.benjamini_hochberg,
                        taxon_identifier=9606):
    network = nx.DiGraph()
    for pathway, name in reactome.get_pathways(taxon_identifier):
        network.add_node(pathway)
        network.nodes[pathway]["pathway"] = name

    for child, parent in reactome.get_pathway_relations():
        if child in network and parent in network:
            network.add_edge(child, parent)

    pathways = {pathway: set() for pathway in network}
    for protein, pathway in reactome.get_pathway_map(taxon_identifier):
        pathways[pathway].add(protein)

    pathways = {
        pathway: pathways[pathway]
        for pathway in pathways if pathways[pathway]
    }
    network.remove_nodes_from(pathway for pathway in network
                              if pathway not in pathways)

    mapped_proteins = set.union(*pathways.values())

    intersection = {
        pathway: len(pathways[pathway].intersection(
            protein_protein_interaction_network.nodes()))
        for pathway in pathways
    }

    p_values = {
        pathway: p
        for pathway, p in correction({
            pathway: test(
                intersection[pathway], len(mapped_proteins),
                len(pathways[pathway]),
                len(
                    mapped_proteins.intersection(
                        protein_protein_interaction_network.nodes())))
            for pathway in network if intersection[pathway]
        }).items()
    }

    for pathway in network:
        network.nodes[pathway]["pathway proteins"] = len(pathways[pathway])
        network.nodes[pathway]["network proteins"] = intersection[pathway]
        network.nodes[pathway]["p-value"] = p_values[pathway]

    return network


def get_pathway_sizes(network):
    return {
        pathway: network.nodes[pathway]["pathway proteins"]
        for pathway in network
    }


def export(network, basename, suffix=""):
    nx.write_graphml_xml(network, "{0}{1}.graphml".format(basename, suffix))
