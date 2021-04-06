import os
import json

import networkx as nx


def export_network(network, path, suffix=""):
    nx.write_graphml_xml(network, "{0}{2}{1}".format(*os.path.splitext(path), suffix))


def export_styles(styles, path, suffix=""):
    styles.write(
        "{0}{2}{1}".format(*os.path.splitext(path), suffix),
        encoding="UTF-8",
        xml_declaration=True,
    )
