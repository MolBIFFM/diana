import networkx as nx
from databases import reactome


def add_pathway_hierarchy(network):
    reactome.add_pathways(network)
    reactome.add_pathway_relations(network)


def export(network, basename, suffix=""):
    nx.write_graphml_xml(network, "{0}{1}.graphml".format(basename, suffix))
