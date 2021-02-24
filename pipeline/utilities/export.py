import os
import json

import networkx as nx


def export_network(path, network, suffix=""):
    if path.endswith(".graphml") or path.endswith(".xml"):
        nx.write_graphml_xml(
            network, "{0}{2}{1}".format(*os.path.splitext(path), suffix)
        )
    elif path.endswith(".cyjs") or path.endswith(".json"):
        with open("{0}{2}{1}".format(*os.path.splitext(path), suffix), "w") as file:
            json.dump(nx.readwrite.json_graph.cytoscape_data(network), file, indent=2)


def export_styles(path, styles):
    if path.endswith(".xml"):
        with open(path, "wb") as file:
            styles.write(file, encoding="UTF-8", xml_declaration=True)
