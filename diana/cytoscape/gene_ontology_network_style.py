"""The Cytoscape style configuration for a Gene Ontology network."""
import math
import os
import xml.etree.ElementTree as ET
from typing import Optional

import networkx as nx
from cytoscape import elements

from networks import gene_ontology_network

COMPONENTS: dict[str, dict[str, dict[str, dict[
    str, bool | int | str | float | dict[str, bool | int | str | float | dict[
        str, bool | int | str |
        float | dict[str, bool | int | str | float]]]]]]] = {
            "edge": {
                "dependency": {
                    "arrowColorMatchesEdge": {
                        "value": True
                    }
                },
                "visualProperty": {
                    "EDGE_BEND": {
                        "default": ""
                    },
                    "EDGE_CURVED": {
                        "default": True
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
                        "default": 10
                    },
                    "EDGE_LABEL_TRANSPARENCY": {
                        "default": 255
                    },
                    "EDGE_LABEL_WIDTH": {
                        "default": 200.0
                    },
                    "EDGE_LINE_TYPE": {
                        "default": "SOLID"
                    },
                    "EDGE_PAINT": {
                        "default": f"#{50:02X}{50:02X}{50:02X}"
                    },
                    "EDGE_SELECTED": {
                        "default": False
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
                        "default": f"#{132:02X}{132:02X}{132:02X}"
                    },
                    "EDGE_TOOLTIP": {
                        "default": ""
                    },
                    "EDGE_TRANSPARENCY": {
                        "default": 255
                    },
                    "EDGE_UNSELECTED_PAINT": {
                        "default": f"#{64:02X}{64:02X}{64:02X}"
                    },
                    "EDGE_VISIBLE": {
                        "default": True
                    },
                    "EDGE_WIDTH": {
                        "default": 2.0
                    },
                },
            },
            "network": {
                "dependency": {
                },
                "visualProperty": {
                    "NETWORK_BACKGROUND_PAINT": {
                        "default": f"#{255:02X}{255:02X}{255:02X}"
                    },
                    "NETWORK_CENTER_X_LOCATION": {
                        "default": 0.0
                    },
                    "NETWORK_CENTER_Y_LOCATION": {
                        "default": 0.0
                    },
                    "NETWORK_CENTER_Z_LOCATION": {
                        "default": 0.0
                    },
                    "NETWORK_DEPTH": {
                        "default": 0.0
                    },
                    "NETWORK_EDGE_SELECTION": {
                        "default": True
                    },
                    "NETWORK_NODE_SELECTION": {
                        "default": True
                    },
                    "NETWORK_SCALE_FACTOR": {
                        "default": 1.0
                    },
                    "NETWORK_SIZE": {
                        "default": 550.0
                    },
                    "NETWORK_TITLE": {
                        "default": ""
                    },
                    "NETWORK_WIDTH": {
                        "default": 550.0
                    },
                },
            },
            "node": {
                "dependency": {
                    "nodeCustomGraphicsSizeSync": {
                        "value": True
                    },
                    "nodeSizeLocked": {
                        "value": True
                    },
                },
                "visualProperty": {
                    "COMPOUND_NODE_PADDING": {
                        "default": 10.0
                    },
                    "NODE_BORDER_PAINT": {
                        "default": f"#{204:02X}{204:02X}{204:02X}"
                    },
                    "NODE_BORDER_STROKE": {
                        "default": "SOLID"
                    },
                    "NODE_BORDER_TRANSPARENCY": {
                        "default": 255
                    },
                    "NODE_BORDER_WIDTH": {
                        "default": 1.0
                    },
                    "NODE_COMPOUND_SHAPE": {
                        "default": "ROUND_RECTANGLE"
                    },
                    "NODE_CUSTOMGRAPHICS_1": {
                        "default": "org.cytoscape.ding.customgraphics."
                                   "NullCustomGraphics,0,[ Remove Graphics ],"
                    },
                    "NODE_CUSTOMGRAPHICS_2": {
                        "default": "org.cytoscape.ding.customgraphics."
                                   "NullCustomGraphics,0,[ Remove Graphics ],"
                    },
                    "NODE_CUSTOMGRAPHICS_3": {
                        "default": "org.cytoscape.ding.customgraphics."
                                   "NullCustomGraphics,0,[ Remove Graphics ],"
                    },
                    "NODE_CUSTOMGRAPHICS_4": {
                        "default": "org.cytoscape.ding.customgraphics."
                                   "NullCustomGraphics,0,[ Remove Graphics ],"
                    },
                    "NODE_CUSTOMGRAPHICS_5": {
                        "default": "org.cytoscape.ding.customgraphics."
                                   "NullCustomGraphics,0,[ Remove Graphics ],"
                    },
                    "NODE_CUSTOMGRAPHICS_6": {
                        "default": "org.cytoscape.ding.customgraphics."
                                   "NullCustomGraphics,0,[ Remove Graphics ],"
                    },
                    "NODE_CUSTOMGRAPHICS_7": {
                        "default": "org.cytoscape.ding.customgraphics."
                                   "NullCustomGraphics,0,[ Remove Graphics ],"
                    },
                    "NODE_CUSTOMGRAPHICS_8": {
                        "default": "org.cytoscape.ding.customgraphics."
                                   "NullCustomGraphics,0,[ Remove Graphics ],"
                    },
                    "NODE_CUSTOMGRAPHICS_9": {
                        "default": "org.cytoscape.ding.customgraphics."
                                   "NullCustomGraphics,0,[ Remove Graphics ],"
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
                        "default": 50.0
                    },
                    "NODE_CUSTOMGRAPHICS_SIZE_2": {
                        "default": 50.0
                    },
                    "NODE_CUSTOMGRAPHICS_SIZE_3": {
                        "default": 50.0
                    },
                    "NODE_CUSTOMGRAPHICS_SIZE_4": {
                        "default": 50.0
                    },
                    "NODE_CUSTOMGRAPHICS_SIZE_5": {
                        "default": 50.0
                    },
                    "NODE_CUSTOMGRAPHICS_SIZE_6": {
                        "default": 50.0
                    },
                    "NODE_CUSTOMGRAPHICS_SIZE_7": {
                        "default": 50.0
                    },
                    "NODE_CUSTOMGRAPHICS_SIZE_8": {
                        "default": 50.0
                    },
                    "NODE_CUSTOMGRAPHICS_SIZE_9": {
                        "default": 50.0
                    },
                    "NODE_CUSTOMPAINT_1": {
                        "default":
                            "DefaultVisualizableVisualProperty("
                            "id=NODE_CUSTOMPAINT_1, name=Node Custom Paint 1)"
                    },
                    "NODE_CUSTOMPAINT_2": {
                        "default":
                            "DefaultVisualizableVisualProperty("
                            "id=NODE_CUSTOMPAINT_2, name=Node Custom Paint 2)"
                    },
                    "NODE_CUSTOMPAINT_3": {
                        "default":
                            "DefaultVisualizableVisualProperty("
                            "id=NODE_CUSTOMPAINT_3, name=Node Custom Paint 3)"
                    },
                    "NODE_CUSTOMPAINT_4": {
                        "default":
                            "DefaultVisualizableVisualProperty("
                            "id=NODE_CUSTOMPAINT_4, name=Node Custom Paint 4)"
                    },
                    "NODE_CUSTOMPAINT_5": {
                        "default":
                            "DefaultVisualizableVisualProperty("
                            "id=NODE_CUSTOMPAINT_5, name=Node Custom Paint 5)"
                    },
                    "NODE_CUSTOMPAINT_6": {
                        "default":
                            "DefaultVisualizableVisualProperty("
                            "id=NODE_CUSTOMPAINT_6, name=Node Custom Paint 6)"
                    },
                    "NODE_CUSTOMPAINT_7": {
                        "default":
                            "DefaultVisualizableVisualProperty("
                            "id=NODE_CUSTOMPAINT_7, name=Node Custom Paint 7)"
                    },
                    "NODE_CUSTOMPAINT_8": {
                        "default":
                            "DefaultVisualizableVisualProperty("
                            "id=NODE_CUSTOMPAINT_8, name=Node Custom Paint 8)"
                    },
                    "NODE_CUSTOMPAINT_9": {
                        "default":
                            "DefaultVisualizableVisualProperty("
                            "id=NODE_CUSTOMPAINT_9, name=Node Custom Paint 9)"
                    },
                    "NODE_DEPTH": {
                        "default": 0.0
                    },
                    "NODE_FILL_COLOR": {
                        "default": f"#{255:02X}{255:02X}{255:02X}",
                    },
                    "NODE_HEIGHT": {
                        "default": 35.0
                    },
                    "NODE_LABEL": {
                        "default": "",
                    },
                    "NODE_LABEL_COLOR": {
                        "default": f"#{0:02X}{0:02X}{0:02X}"
                    },
                    "NODE_LABEL_FONT_FACE": {
                        "default": "SansSerif.plain,plain,12"
                    },
                    "NODE_LABEL_FONT_SIZE": {
                        "default": 12
                    },
                    "NODE_LABEL_POSITION": {
                        "default": "C,C,c,0.00,0.00"
                    },
                    "NODE_LABEL_TRANSPARENCY": {
                        "default": 255
                    },
                    "NODE_LABEL_WIDTH": {
                        "default": 200.0
                    },
                    "NODE_NESTED_NETWORK_IMAGE_VISIBLE": {
                        "default": True
                    },
                    "NODE_PAINT": {
                        "default": f"#{30:02X}{144:02X}{255:02X}"
                    },
                    "NODE_SELECTED": {
                        "default": False
                    },
                    "NODE_SELECTED_PAINT": {
                        "default": f"#{255:02X}{255:02X}{0:02X}"
                    },
                    "NODE_SHAPE": {
                        "default": "ROUND_RECTANGLE"
                    },
                    "NODE_SIZE": {
                        "default": 35.0
                    },
                    "NODE_TOOLTIP": {
                        "default": ""
                    },
                    "NODE_TRANSPARENCY": {
                        "default": 255
                    },
                    "NODE_WIDTH": {
                        "default": 75.0
                    },
                    "NODE_X_LOCATION": {
                        "default": 0.0
                    },
                    "NODE_Y_LOCATION": {
                        "default": 0.0
                    },
                    "NODE_Z_LOCATION": {
                        "default": 0.0
                    },
                },
            },
        }


def get_styles(network: nx.Graph) -> ET.ElementTree:
    """
    Returns the Cytoscape styles for a Gene Ontology network.

    Args:
        network: The Gene Ontology network.

    Returns:
        The Cytoscape style for the Gene Ontology network.
    """
    # Compile the element tree for style specification of network components.
    styles = ET.ElementTree(
        ET.Element("vizmap", attrib={
            "id": "VizMap",
            "documentVersion": "3.0"
        }))

    visual_style_sub_element = ET.SubElement(styles.getroot(),
                                             "visualStyle",
                                             attrib={"name": "Gene Ontology"})

    visual_properties: dict[str, dict[str, ET.Element]] = {}
    for component, properties in COMPONENTS.items():
        visual_properties[component] = {}
        component_sub_element = ET.SubElement(visual_style_sub_element,
                                              component)
        for name, dependency in properties["dependency"].items():
            if isinstance(dependency["value"], bool):
                elements.add_dependency(component_sub_element, name,
                                        dependency["value"])

        for name, visual_property in properties["visualProperty"].items():
            if isinstance(visual_property["default"], (bool, int, str, float)):
                visual_properties[component][
                    name] = elements.add_visual_property(
                        component_sub_element, name, visual_property["default"])

    # Assign the Gene Ontology term identifier as node label.
    elements.add_passthrough_mapping(visual_properties["node"]["NODE_LABEL"],
                                     "name", "string")

    # Assign the Gene Ontology term description as node tooltip.
    elements.add_passthrough_mapping(visual_properties["node"]["NODE_TOOLTIP"],
                                     "term", "string")

    # Assign node shape according to the Gene Ontology namespace.
    elements.add_discrete_mapping(
        visual_properties["node"]["NODE_SHAPE"], "namespace", "string", {
            "cellular component": "RECTANGLE",
            "biological process": "TRIANGLE",
            "molecular function": "ELLIPSE",
        })

    # Color nodes red as a linear function of p-value.
    elements.add_continuous_mapping(
        visual_properties["node"]["NODE_FILL_COLOR"], "p-value", "float", {
            0.0: (f"#{255:02X}{0:02X}{0:02X}", f"#{255:02X}{0:02X}{0:02X}",
                  f"#{255:02X}{0:02X}{0:02X}"),
            1.0: (f"#{255:02X}{255:02X}{255:02X}",
                  f"#{255:02X}{255:02X}{255:02X}",
                  f"#{255:02X}{255:02X}{255:02X}")
        })

    # Scale nodes as a function of the number of associated proteins.
    max_number_proteins = max(
        gene_ontology_network.get_term_sizes(network).values())

    if isinstance(COMPONENTS["node"]["visualProperty"]["NODE_SIZE"]["default"],
                  float):
        elements.add_continuous_mapping(
            visual_properties["node"]["NODE_SIZE"],
            "number of associated proteins", "integer", {
                0: (COMPONENTS["node"]["visualProperty"]["NODE_SIZE"]
                    ["default"], COMPONENTS["node"]["visualProperty"]
                    ["NODE_SIZE"]["default"],
                    COMPONENTS["node"]["visualProperty"]["NODE_SIZE"]["default"]
                   ),
                max_number_proteins:
                    (COMPONENTS["node"]["visualProperty"]["NODE_SIZE"]
                     ["default"] + math.sqrt(max_number_proteins),) * 3
            })

    # Indent the representation of the element tree.
    ET.indent(styles)

    # Return the Cytoscape style for the Gene Ontology network.
    return styles


def export(styles: ET.ElementTree, basename: str) -> Optional[str]:
    """
    Exports the Cytoscape style for a Gene Ontology network if possible without
    overwriting an existing file.

    Args:
        styles: The Cytoscape styles.
        basename: The base file name.

    Returns:
        The file name the Cytoscape styles were exported to if there is no
        naming conflict.
    """
    # Avoid overwriting an existing file.
    if os.path.isfile(f"{basename}.xml"):
        return None

    # Export the element tree of the Cytoscape style.
    styles.write(
        f"{basename}.xml",
        encoding="utf-8",
        xml_declaration=True,
    )

    # Return the file name of the Cytoscape style.
    return f"{basename}.xml"
