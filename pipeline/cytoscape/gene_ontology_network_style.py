"""The Cytoscape style configuration for a Gene Ontology network."""

import xml.etree.ElementTree as ET

import networkx as nx
from networks import gene_ontology_network

COMPONENTS = {
    "edge": {
        "dependency": {
            "arrowColorMatchesEdge": {
                "value": "false"
            }
        },
        "visualProperty": {
            "EDGE_BEND": {
                "default": ""
            },
            "EDGE_CURVED": {
                "default": "true"
            },
            "EDGE_LABEL": {
                "default": ""
            },
            "EDGE_LABEL_COLOR": {
                "default": f"#{0:02X}{0:02X}{0:02X}"
            },
            "EDGE_LABEL_FONT_FACE": {
                "default": "Dialog.plain,plain,10"
            },
            "EDGE_LABEL_FONT_SIZE": {
                "default": "10"
            },
            "EDGE_LABEL_TRANSPARENCY": {
                "default": "255"
            },
            "EDGE_LABEL_WIDTH": {
                "default": "200.0"
            },
            "EDGE_LINE_TYPE": {
                "default": "SOLID"
            },
            "EDGE_PAINT": {
                "default": f"#{50:02X}{50:02X}{50:02X}"
            },
            "EDGE_SELECTED": {
                "default": "false"
            },
            "EDGE_SELECTED_PAINT": {
                "default": f"#{255:02X}{0:02X}{0:02X}"
            },
            "EDGE_SOURCE_ARROW_SELECTED_PAINT": {
                "default": f"#{255:02X}{255:02X}{0:02X}"
            },
            "EDGE_SOURCE_ARROW_SHAPE": {
                "default": "NONE"
            },
            "EDGE_SOURCE_ARROW_UNSELECTED_PAINT": {
                "default": f"#{0:02X}{0:02X}{0:02X}"
            },
            "EDGE_STROKE_SELECTED_PAINT": {
                "default": f"#{255:02X}{0:02X}{0:02X}"
            },
            "EDGE_STROKE_UNSELECTED_PAINT": {
                "default": f"#{132:02X}{132:02X}{132:02X}"
            },
            "EDGE_TARGET_ARROW_SELECTED_PAINT": {
                "default": f"#{255:02X}{255:02X}{0:02X}"
            },
            "EDGE_TARGET_ARROW_SHAPE": {
                "default": "DELTA"
            },
            "EDGE_TARGET_ARROW_UNSELECTED_PAINT": {
                "default": f"#{0:02X}{0:02X}{0:02X}"
            },
            "EDGE_TOOLTIP": {
                "default": ""
            },
            "EDGE_TRANSPARENCY": {
                "default": "255"
            },
            "EDGE_UNSELECTED_PAINT": {
                "default": f"#{64:02X}{64:02X}{64:02X}"
            },
            "EDGE_VISIBLE": {
                "default": "true"
            },
            "EDGE_WIDTH": {
                "default": "2.0"
            },
        },
    },
    "network": {
        "dependency": {},
        "visualProperty": {
            "NETWORK_BACKGROUND_PAINT": {
                "default": f"#{255:02X}{255:02X}{255:02X}"
            },
            "NETWORK_CENTER_X_LOCATION": {
                "default": "0.0"
            },
            "NETWORK_CENTER_Y_LOCATION": {
                "default": "0.0"
            },
            "NETWORK_CENTER_Z_LOCATION": {
                "default": "0.0"
            },
            "NETWORK_DEPTH": {
                "default": "0.0"
            },
            "NETWORK_EDGE_SELECTION": {
                "default": "true"
            },
            "NETWORK_NODE_SELECTION": {
                "default": "true"
            },
            "NETWORK_SCALE_FACTOR": {
                "default": "1.0"
            },
            "NETWORK_SIZE": {
                "default": "550.0"
            },
            "NETWORK_TITLE": {
                "default": ""
            },
            "NETWORK_WIDTH": {
                "default": "550.0"
            },
        },
    },
    "node": {
        "dependency": {
            "nodeCustomGraphicsSizeSync": {
                "value": "true"
            },
            "nodeSizeLocked": {
                "value": "true"
            },
        },
        "visualProperty": {
            "COMPOUND_NODE_PADDING": {
                "default": "10.0"
            },
            "NODE_BORDER_PAINT": {
                "default": f"#{204:02X}{204:02X}{204:02X}"
            },
            "NODE_BORDER_STROKE": {
                "default": "SOLID"
            },
            "NODE_BORDER_TRANSPARENCY": {
                "default": "255"
            },
            "NODE_BORDER_WIDTH": {
                "default": "0.0"
            },
            "NODE_COMPOUND_SHAPE": {
                "default": "ROUND_RECTANGLE"
            },
            "NODE_CUSTOMGRAPHICS_1": {
                "default":
                    "org.cytoscape.ding.customgraphics.NullCustomGraphics,0,"
                    "[ Remove Graphics ],"
            },
            "NODE_CUSTOMGRAPHICS_2": {
                "default":
                    "org.cytoscape.ding.customgraphics.NullCustomGraphics,0,"
                    "[ Remove Graphics ],"
            },
            "NODE_CUSTOMGRAPHICS_3": {
                "default":
                    "org.cytoscape.ding.customgraphics.NullCustomGraphics,0,"
                    "[ Remove Graphics ],"
            },
            "NODE_CUSTOMGRAPHICS_4": {
                "default":
                    "org.cytoscape.ding.customgraphics.NullCustomGraphics,0,"
                    "[ Remove Graphics ],"
            },
            "NODE_CUSTOMGRAPHICS_5": {
                "default":
                    "org.cytoscape.ding.customgraphics.NullCustomGraphics,0,"
                    "[ Remove Graphics ],"
            },
            "NODE_CUSTOMGRAPHICS_6": {
                "default":
                    "org.cytoscape.ding.customgraphics.NullCustomGraphics,0,"
                    "[ Remove Graphics ],"
            },
            "NODE_CUSTOMGRAPHICS_7": {
                "default":
                    "org.cytoscape.ding.customgraphics.NullCustomGraphics,0,"
                    "[ Remove Graphics ],"
            },
            "NODE_CUSTOMGRAPHICS_8": {
                "default":
                    "org.cytoscape.ding.customgraphics.NullCustomGraphics,0,"
                    "[ Remove Graphics ],"
            },
            "NODE_CUSTOMGRAPHICS_9": {
                "default":
                    "org.cytoscape.ding.customgraphics.NullCustomGraphics,0,"
                    "[ Remove Graphics ],"
            },
            "NODE_CUSTOMGRAPHICS_POSITION_1": {
                "default": "C,C,c,0.00,0.00"
            },
            "NODE_CUSTOMGRAPHICS_POSITION_2": {
                "default": "C,C,c,0.00,0.00"
            },
            "NODE_CUSTOMGRAPHICS_POSITION_3": {
                "default": "C,C,c,0.00,0.00"
            },
            "NODE_CUSTOMGRAPHICS_POSITION_4": {
                "default": "C,C,c,0.00,0.00"
            },
            "NODE_CUSTOMGRAPHICS_POSITION_5": {
                "default": "C,C,c,0.00,0.00"
            },
            "NODE_CUSTOMGRAPHICS_POSITION_6": {
                "default": "C,C,c,0.00,0.00"
            },
            "NODE_CUSTOMGRAPHICS_POSITION_7": {
                "default": "C,C,c,0.00,0.00"
            },
            "NODE_CUSTOMGRAPHICS_POSITION_8": {
                "default": "C,C,c,0.00,0.00"
            },
            "NODE_CUSTOMGRAPHICS_POSITION_9": {
                "default": "C,C,c,0.00,0.00"
            },
            "NODE_CUSTOMGRAPHICS_SIZE_1": {
                "default": "50.0"
            },
            "NODE_CUSTOMGRAPHICS_SIZE_2": {
                "default": "50.0"
            },
            "NODE_CUSTOMGRAPHICS_SIZE_3": {
                "default": "50.0"
            },
            "NODE_CUSTOMGRAPHICS_SIZE_4": {
                "default": "50.0"
            },
            "NODE_CUSTOMGRAPHICS_SIZE_5": {
                "default": "50.0"
            },
            "NODE_CUSTOMGRAPHICS_SIZE_6": {
                "default": "50.0"
            },
            "NODE_CUSTOMGRAPHICS_SIZE_7": {
                "default": "50.0"
            },
            "NODE_CUSTOMGRAPHICS_SIZE_8": {
                "default": "50.0"
            },
            "NODE_CUSTOMGRAPHICS_SIZE_9": {
                "default": "50.0"
            },
            "NODE_CUSTOMPAINT_1": {
                "default":
                    "DefaultVisualizableVisualProperty(id=NODE_CUSTOMPAINT_1, "
                    "name=Node Custom Paint 1)"
            },
            "NODE_CUSTOMPAINT_2": {
                "default":
                    "DefaultVisualizableVisualProperty(id=NODE_CUSTOMPAINT_2, "
                    "name=Node Custom Paint 2)"
            },
            "NODE_CUSTOMPAINT_3": {
                "default":
                    "DefaultVisualizableVisualProperty(id=NODE_CUSTOMPAINT_3, "
                    "name=Node Custom Paint 3)"
            },
            "NODE_CUSTOMPAINT_4": {
                "default":
                    "DefaultVisualizableVisualProperty(id=NODE_CUSTOMPAINT_4, "
                    "name=Node Custom Paint 4)"
            },
            "NODE_CUSTOMPAINT_5": {
                "default":
                    "DefaultVisualizableVisualProperty(id=NODE_CUSTOMPAINT_5, "
                    "name=Node Custom Paint 5)"
            },
            "NODE_CUSTOMPAINT_6": {
                "default":
                    "DefaultVisualizableVisualProperty(id=NODE_CUSTOMPAINT_6, "
                    "name=Node Custom Paint 6)"
            },
            "NODE_CUSTOMPAINT_7": {
                "default":
                    "DefaultVisualizableVisualProperty(id=NODE_CUSTOMPAINT_7, "
                    "name=Node Custom Paint 7)"
            },
            "NODE_CUSTOMPAINT_8": {
                "default":
                    "DefaultVisualizableVisualProperty(id=NODE_CUSTOMPAINT_8, "
                    "name=Node Custom Paint 8)"
            },
            "NODE_CUSTOMPAINT_9": {
                "default":
                    "DefaultVisualizableVisualProperty(id=NODE_CUSTOMPAINT_9, "
                    "name=Node Custom Paint 9)"
            },
            "NODE_DEPTH": {
                "default": "0.0"
            },
            "NODE_FILL_COLOR": {
                "default": f"#{137:02X}{208:02X}{255:02X}",
                "continuousMapping": {
                    "attributeName": "p-value",
                    "attributeType": "float",
                    "continuousMappingPoint": {
                        "0.0": {
                            "equalValue": f"#{255:02X}{0:02X}{0:02X}",
                            "greaterValue": f"#{255:02X}{0:02X}{0:02X}",
                            "lesserValue": f"#{255:02X}{0:02X}{0:02X}"
                        },
                        "1.0": {
                            "equalValue": f"#{255:02X}{255:02X}{255:02X}",
                            "greaterValue": f"#{255:02X}{255:02X}{255:02X}",
                            "lesserValue": f"#{255:02X}{255:02X}{255:02X}"
                        }
                    }
                }
            },
            "NODE_HEIGHT": {
                "default": "35.0"
            },
            "NODE_LABEL": {
                "default": "",
                "passthroughMapping": {
                    "attributeName": "name",
                    "attributeType": "string",
                },
            },
            "NODE_LABEL_COLOR": {
                "default": f"#{0:02X}{0:02X}{0:02X}"
            },
            "NODE_LABEL_FONT_FACE": {
                "default": "SansSerif.plain,plain,12"
            },
            "NODE_LABEL_FONT_SIZE": {
                "default": "12"
            },
            "NODE_LABEL_POSITION": {
                "default": "C,C,c,0.00,0.00"
            },
            "NODE_LABEL_TRANSPARENCY": {
                "default": "255"
            },
            "NODE_LABEL_WIDTH": {
                "default": "200.0"
            },
            "NODE_NESTED_NETWORK_IMAGE_VISIBLE": {
                "default": "true"
            },
            "NODE_PAINT": {
                "default": f"#{30:02X}{144:02X}{255:02X}"
            },
            "NODE_SELECTED": {
                "default": "false"
            },
            "NODE_SELECTED_PAINT": {
                "default": f"#{255:02X}{255:02X}{0:02X}"
            },
            "NODE_SHAPE": {
                "default": "ROUND_RECTANGLE",
                "discreteMapping": {
                    "attributeName": "namespace",
                    "attributeType": "string",
                    "discreteMappingEntry": {
                        "cellular_component": "RECTANGLE",
                        "biological_process": "TRIANGLE",
                        "molecular_function": "ELLIPSE",
                    },
                },
            },
            "NODE_SIZE": {
                "default": "35.0",
                "continuousMapping": {
                    "attributeName": "proteins",
                    "attributeType": "float",
                    "continuousMappingPoint": {
                        "1.0": {
                            "equalValue": "36.0",
                            "greaterValue": "36.0",
                            "lesserValue": "36.0"
                        },
                        "{max}": {
                            "equalValue": "{size}",
                            "greaterValue": "{size}",
                            "lesserValue": "{size}"
                        },
                    }
                }
            },
            "NODE_TOOLTIP": {
                "default": "",
                "passthroughMapping": {
                    "attributeName": "term",
                    "attributeType": "string",
                },
            },
            "NODE_TRANSPARENCY": {
                "default": "255"
            },
            "NODE_WIDTH": {
                "default": "75.0"
            },
            "NODE_X_LOCATION": {
                "default": "0.0"
            },
            "NODE_Y_LOCATION": {
                "default": "0.0"
            },
            "NODE_Z_LOCATION": {
                "default": "0.0"
            },
        },
    },
}


def get_gene_ontology_network_style(network: nx.Graph) -> ET.ElementTree:
    """
    Returns the Cytoscape styles for a Gene Ontology network.

    Args:
        network: The Gene Ontology network.

    Returns:
        The Cytoscape style for the Gene Ontology network.
    """
    style = ET.ElementTree(
        ET.Element("vizmap", attrib={
            "id": "VizMap",
            "documentVersion": "3.0"
        }))

    visual_style_sub_element = ET.SubElement(style.getroot(),
                                             "visualStyle",
                                             attrib={"name": "ontology"})

    for component in COMPONENTS:
        component_sub_element = ET.SubElement(visual_style_sub_element,
                                              component)

        for name, dependency in COMPONENTS[component]["dependency"].items():
            ET.SubElement(
                component_sub_element,
                "dependency",
                attrib={
                    "name": name,
                    "value": dependency["value"]
                },
            )

        for name, visual_property in COMPONENTS[component][
                "visualProperty"].items():
            visual_property_sub_element = ET.SubElement(
                component_sub_element,
                "visualProperty",
                attrib={
                    "name": name,
                    "default": visual_property["default"]
                },
            )

            if visual_property.get("continuousMapping"):
                continuous_mapping_sub_element = ET.SubElement(
                    visual_property_sub_element,
                    "continuousMapping",
                    attrib={
                        "attributeName":
                            visual_property["continuousMapping"]
                            ["attributeName"],
                        "attributeType":
                            visual_property["continuousMapping"]
                            ["attributeType"]
                    })

                for key, values in visual_property["continuousMapping"][
                        "continuousMappingPoint"].items():
                    ET.SubElement(
                        continuous_mapping_sub_element,
                        "continuousMappingPoint",
                        attrib={
                            "attrValue":
                                key.format(max=float(
                                    max(
                                        gene_ontology_network.get_term_sizes(
                                            network).values()))),
                            "equalValue":
                                values["equalValue"].format(size=35.0 + float(
                                    max(
                                        gene_ontology_network.get_term_sizes(
                                            network).values()))),
                            "greaterValue":
                                values["greaterValue"].format(size=35.0 + float(
                                    max(
                                        gene_ontology_network.get_term_sizes(
                                            network).values()))),
                            "lesserValue":
                                values["lesserValue"].format(size=35.0 + float(
                                    max(
                                        gene_ontology_network.get_term_sizes(
                                            network).values())))
                        },
                    )

            elif visual_property.get("discreteMapping"):
                discrete_mapping_sub_element = ET.SubElement(
                    visual_property_sub_element,
                    "discreteMapping",
                    attrib={
                        "attributeName":
                            visual_property["discreteMapping"]["attributeName"],
                        "attributeType":
                            visual_property["discreteMapping"]["attributeType"],
                    },
                )

                for key, value in visual_property["discreteMapping"][
                        "discreteMappingEntry"].items():
                    ET.SubElement(
                        discrete_mapping_sub_element,
                        "discreteMappingEntry",
                        attrib={
                            "attributeValue": key,
                            "value": value,
                        },
                    )

            elif visual_property.get("passthroughMapping"):
                ET.SubElement(
                    visual_property_sub_element,
                    "passthroughMapping",
                    attrib={
                        "attributeName":
                            visual_property["passthroughMapping"]
                            ["attributeName"],
                        "attributeType":
                            visual_property["passthroughMapping"]
                            ["attributeType"],
                    },
                )

    ET.indent(style)
    return style


def export(styles: ET.ElementTree, basename: str) -> None:
    """
    Exports the Cytoscape style.

    Args:
        styles: The Cytoscape style.
        basename: The base file name.
    """
    styles.write(
        f"{basename}.xml",
        encoding="utf-8",
        xml_declaration=True,
    )
