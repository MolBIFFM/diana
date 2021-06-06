import networkx as nx


def export_network(network, basename, suffix=""):
    nx.write_graphml_xml(network, "{0}{1}.graphml".format(basename, suffix))


def export_styles(styles, basename, suffix=""):
    styles.write(
        "{0}{1}.xml".format(basename, suffix),
        encoding="UTF-8",
        xml_declaration=True,
    )
