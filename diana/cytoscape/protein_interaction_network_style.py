"""
The Cytoscape style configuration for a protein-protein interaction
network.
"""

import json
import statistics
import xml.etree.ElementTree as ET
from typing import Callable, Collection, Iterable

import networkx as nx
from networks import protein_interaction_network

COMPONENTS: list[dict[str, dict[str, dict[str, dict[
    str, str | dict[str, str | dict[str, str | dict[str, str]]]]]]]] = [
        {
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
                        "default": "NONE"
                    },
                    "EDGE_TARGET_ARROW_UNSELECTED_PAINT": {
                        "default": f"#{0:02X}{0:02X}{0:02X}"
                    },
                    "EDGE_TOOLTIP": {
                        "default": ""
                    },
                    "EDGE_TRANSPARENCY": {
                        "default": "255",
                        "continuousMapping": {
                            "attributeName": "score",
                            "attributeType": "float",
                            "continuousMappingPoint": {
                                "0": {
                                    "equalValue": "0",
                                    "greaterValue": "0",
                                    "lesserValue": "0"
                                },
                                "{max:.1f}": {
                                    "equalValue": "255",
                                    "greaterValue": "255",
                                    "lesserValue": "255"
                                },
                            }
                        }
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
                        "value": "false"
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
                        "default": "0.0"
                    },
                    "NODE_FILL_COLOR": {
                        "default": f"#{128:02X}{128:02X}{128:02X}",
                        "discreteMapping": {
                            "attributeName": "measurement {time}",
                            "attributeType": "string",
                            "discreteMappingEntry": {
                                "up": f"#{255:02X}{0:02X}{0:02X}",
                                "mid up": f"#{255:02X}{102:02X}{102:02X}",
                                "mid": f"#{128:02X}{128:02X}{128:02X}",
                                "mid down": f"#{102:02X}{102:02X}{255:02X}",
                                "down": f"#{0:02X}{0:02X}{255:02X}",
                            },
                        },
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
                            "attributeName":
                                "post-translational modification {time}",
                            "attributeType":
                                "string",
                            "discreteMappingEntry": {
                                "{modifications[0]}": "RECTANGLE",
                            },
                        },
                    },
                    "NODE_SIZE": {
                        "default": "35.0"
                    },
                    "NODE_TOOLTIP": {
                        "default": "",
                        "passthroughMapping": {
                            "attributeName": "protein",
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
        },
        {
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
                        "default": "NONE"
                    },
                    "EDGE_TARGET_ARROW_UNSELECTED_PAINT": {
                        "default": f"#{0:02X}{0:02X}{0:02X}"
                    },
                    "EDGE_TOOLTIP": {
                        "default": ""
                    },
                    "EDGE_TRANSPARENCY": {
                        "default": "255",
                        "continuousMapping": {
                            "attributeName": "score",
                            "attributeType": "float",
                            "continuousMappingPoint": {
                                "0.0": {
                                    "equalValue": "0",
                                    "greaterValue": "0",
                                    "lesserValue": "0"
                                },
                                "{max:.1f}": {
                                    "equalValue": "255",
                                    "greaterValue": "255",
                                    "lesserValue": "255"
                                },
                            }
                        }
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
                        "value": "false"
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
                        "default": "0.0"
                    },
                    "NODE_FILL_COLOR": {
                        "default": f"#{128:02X}{128:02X}{128:02X}",
                        "discreteMapping": {
                            "attributeName": "measurement {time}",
                            "attributeType": "string",
                            "discreteMappingEntry": {
                                "up":
                                    f"#{255:02X}{0:02X}{0:02X}",
                                "mid up":
                                    f"#{255:02X}{102:02X}{102:02X}",
                                "mid":
                                    f"#{128:02X}{128:02X}{128:02X}",
                                "mid down":
                                    f"#{102:02X}{102:02X}{255:02X}",
                                "down":
                                    f"#{0:02X}{0:02X}{255:02X}",
                                "{modifications[0]} up "
                                "{modifications[1]} down":
                                    f"#{0:02X}{255:02X}{0:02X}",
                                "{modifications[0]} up "
                                "{modifications[1]} mid down":
                                    f"#{0:02X}{255:02X}{0:02X}",
                                "{modifications[0]} mid up "
                                "{modifications[1]} down":
                                    f"#{0:02X}{255:02X}{0:02X}",
                                "{modifications[0]} mid up "
                                "{modifications[1]} mid down":
                                    f"#{0:02X}{255:02X}{0:02X}",
                                "{modifications[0]} down "
                                "{modifications[1]} up":
                                    f"#{255:02X}{255:02X}{0:02X}",
                                "{modifications[0]} down "
                                "{modifications[1]} mid up":
                                    f"#{255:02X}{255:02X}{0:02X}",
                                "{modifications[0]} mid down "
                                "{modifications[1]} up":
                                    f"#{255:02X}{255:02X}{0:02X}",
                                "{modifications[0]} mid down "
                                "{modifications[1]} mid up":
                                    f"#{255:02X}{255:02X}{0:02X}",
                                "{modifications[0]} up "
                                "{modifications[1]} mid":
                                    f"#{255:02X}{0:02X}{0:02X}",
                                "{modifications[0]} mid up "
                                "{modifications[1]} mid":
                                    f"#{255:02X}{102:02X}{102:02X}",
                                "{modifications[0]} mid down "
                                "{modifications[1]} mid":
                                    f"#{102:02X}{102:02X}{255:02X}",
                                "{modifications[0]} down "
                                "{modifications[1]} mid":
                                    f"#{0:02X}{0:02X}{255:02X}",
                                "{modifications[0]} mid "
                                "{modifications[1]} up":
                                    f"#{255:02X}{0:02X}{0:02X}",
                                "{modifications[0]} mid "
                                "{modifications[1]} mid up":
                                    f"#{255:02X}{102:02X}{102:02X}",
                                "{modifications[0]} mid "
                                "{modifications[1]} mid down":
                                    f"#{102:02X}{102:02X}{255:02X}",
                                "{modifications[0]} mid "
                                "{modifications[1]} down":
                                    f"#{0:02X}{0:02X}{255:02X}",
                            },
                        },
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
                            "attributeName":
                                "post-translational modification {time}",
                            "attributeType":
                                "string",
                            "discreteMappingEntry": {
                                "{modifications[0]}":
                                    "RECTANGLE",
                                "{modifications[1]}":
                                    "TRIANGLE",
                                "{modifications[0]} {modifications[1]}":
                                    "ELLIPSE",
                            },
                        },
                    },
                    "NODE_SIZE": {
                        "default": "35.0"
                    },
                    "NODE_TOOLTIP": {
                        "default": "",
                        "passthroughMapping": {
                            "attributeName": "protein",
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
        },
    ]


def get_bar_chart(
    time: int,
    modification: str,
    sites: int,
    cy_range: tuple[float, float] = (-2.0, 2.0)) -> str:
    """
    Returns a bar chart specification for Cytoscape styles.

    Args:
        time: time of measurement associated with the measurement
        modification: modification associated with the measurement
        cy_range: range of binary logarithms of measurements covered by the bar
            chart.

    Returns:
       bar chart specification
    """
    bar_chart = json.dumps({
        "cy_range":
            cy_range,
        "cy_showRangeAxis":
            True,
        "cy_type":
            "UP_DOWN",
        "cy_autoRange":
            False,
        "cy_colorScheme":
            "BLUE_RED",
        "cy_showRangeZeroBaseline":
            True,
        "cy_colors": ["#FF0000", "#0000FF"],
        "cy_dataColumns": [
            f"{time} {modification} {site}" for site in range(1, sites + 1)
        ],
    })
    return f"org.cytoscape.BarChart: {bar_chart}"


def get_style(
    network: nx.Graph,
    bar_chart_range: tuple[float, float] = (-1.0, 1.0),
    convert_measurement: Callable[
        [float, Collection[float]],
        float] = lambda measurement, measurements: measurement,
    site_combination: Callable[[Iterable[float]],
                               float] = lambda sites: max(sites, key=abs),
    replicate_combination: Callable[[Iterable[float]], float] = statistics.mean,
    confidence_score_combination: Callable[[dict[str, float]], float] = lambda
    confidence_scores: float(bool(confidence_scores.values()))
) -> ET.ElementTree:
    """
    Returns the Cytoscape styles for a protein-protein interaction network.

    Args:
        network: The protein-protein interaction network.
        bar_chart_range: The range of binary logarithms of measurements covered
            by the bar chart.
        convert_measurement: The function to transform binary logarithms of
            measurements.
        site_combination: The function to derive a protein-specific measurement
            from its site-specific measurements.
        replicate_combination: The function to derive a site-specific
            measurement from its replicates.
        confidence_score_combination: The function to derive an edge-specific
            confidence score from confidence scores reported in different
            databases.

    Returns:
       The Cytoscape styles for the protein-protein interaction network.
    """
    styles = ET.ElementTree(
        ET.Element("vizmap", attrib={
            "id": "VizMap",
            "documentVersion": "3.0"
        }))

    max_edge_score = confidence_score_combination({
        database: 1.0
        for database in protein_interaction_network.get_databases(network)
    })

    for time in protein_interaction_network.get_times(network):
        modifications = protein_interaction_network.get_modifications(
            network, time)
        visual_style_sub_element = ET.SubElement(styles.getroot(),
                                                 "visualStyle",
                                                 attrib={"name": str(time)})

        for component in COMPONENTS[len(modifications) - 1]:
            component_sub_element = ET.SubElement(visual_style_sub_element,
                                                  component)

            for name, dependency in COMPONENTS[
                    len(modifications) - 1][component]["dependency"].items():
                if isinstance(dependency["value"], str):
                    ET.SubElement(
                        component_sub_element,
                        "dependency",
                        attrib={
                            "name": name,
                            "value": dependency["value"]
                        },
                    )

            for name, visual_property in COMPONENTS[
                    len(modifications) -
                    1][component]["visualProperty"].items():
                if isinstance(visual_property["default"], str):
                    visual_property_sub_element = ET.SubElement(
                        component_sub_element,
                        "visualProperty",
                        attrib={
                            "name": name,
                            "default": visual_property["default"]
                        },
                    )

                if visual_property.get("continuousMapping"):
                    if isinstance(visual_property["continuousMapping"], dict):
                        if isinstance(
                                visual_property["continuousMapping"]
                            ["attributeName"], str) and isinstance(
                                visual_property["continuousMapping"]
                                ["attributeType"], str):
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

                        if isinstance(
                                visual_property["continuousMapping"]
                            ["continuousMappingPoint"], dict):
                            for key, values in visual_property[
                                    "continuousMapping"][
                                        "continuousMappingPoint"].items():
                                if isinstance(values, dict):
                                    if name == "EDGE_TRANSPARENCY":
                                        ET.SubElement(
                                            continuous_mapping_sub_element,
                                            "continuousMappingPoint",
                                            attrib={
                                                "attrValue":
                                                    key.format(
                                                        max=max_edge_score),
                                                "equalValue":
                                                    values["equalValue"],
                                                "greaterValue":
                                                    values["greaterValue"],
                                                "lesserValue":
                                                    values["lesserValue"]
                                            },
                                        )
                                    else:
                                        ET.SubElement(
                                            continuous_mapping_sub_element,
                                            "continuousMappingPoint",
                                            attrib={
                                                "attrValue":
                                                    key,
                                                "equalValue":
                                                    values["equalValue"],
                                                "greaterValue":
                                                    values["greaterValue"],
                                                "lesserValue":
                                                    values["lesserValue"]
                                            },
                                        )

                elif visual_property.get("discreteMapping"):
                    if isinstance(visual_property["discreteMapping"], dict):
                        if isinstance(
                                visual_property["discreteMapping"]
                            ["attributeName"], str) and isinstance(
                                visual_property["discreteMapping"]
                                ["attributeType"], str):
                            if name in ("NODE_FILL_COLOR", "NODE_SHAPE"):
                                discrete_mapping_sub_element = ET.SubElement(
                                    visual_property_sub_element,
                                    "discreteMapping",
                                    attrib={
                                        "attributeName":
                                            visual_property["discreteMapping"]
                                            ["attributeName"].format(time=time),
                                        "attributeType":
                                            visual_property["discreteMapping"]
                                            ["attributeType"],
                                    },
                                )
                            else:
                                discrete_mapping_sub_element = ET.SubElement(
                                    visual_property_sub_element,
                                    "discreteMapping",
                                    attrib={
                                        "attributeName":
                                            visual_property["discreteMapping"]
                                            ["attributeName"],
                                        "attributeType":
                                            visual_property["discreteMapping"]
                                            ["attributeType"],
                                    },
                                )

                        if isinstance(
                                visual_property["discreteMapping"]
                            ["discreteMappingEntry"], dict):
                            for key, value in visual_property[
                                    "discreteMapping"][
                                        "discreteMappingEntry"].items():
                                if isinstance(key, str) and isinstance(
                                        value, str):
                                    if name in ("NODE_FILL_COLOR",
                                                "NODE_SHAPE"):
                                        ET.SubElement(
                                            discrete_mapping_sub_element,
                                            "discreteMappingEntry",
                                            attrib={
                                                "attributeValue":
                                                    key.format(modifications=
                                                               modifications),
                                                "value":
                                                    value,
                                            },
                                        )
                                    else:
                                        ET.SubElement(
                                            discrete_mapping_sub_element,
                                            "discreteMappingEntry",
                                            attrib={
                                                "attributeValue": key,
                                                "value": value,
                                            },
                                        )

                elif visual_property.get("passthroughMapping"):
                    if isinstance(visual_property["passthroughMapping"], dict):
                        if isinstance(
                                visual_property["passthroughMapping"]
                            ["attributeName"], str) and isinstance(
                                visual_property["passthroughMapping"]
                                ["attributeType"], str):
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

        if node_sub_elements := visual_style_sub_element.find("node"):
            visual_property_sub_elements = node_sub_elements.findall(
                "visualProperty")

            for i, modification in enumerate(modifications[:2]):
                for visual_property_sub_element in visual_property_sub_elements:
                    if visual_property_sub_element.get(
                            "name") == f"NODE_CUSTOMGRAPHICS_{i+1}":
                        visual_property_sub_element.set(
                            "default",
                            get_bar_chart(
                                time,
                                modification,
                                protein_interaction_network.get_sites(
                                    network, time, modification),
                                cy_range=(convert_measurement(
                                    bar_chart_range[0],
                                    protein_interaction_network.
                                    get_measurements(network, time,
                                                     modification,
                                                     site_combination,
                                                     replicate_combination)),
                                          convert_measurement(
                                              bar_chart_range[1],
                                              protein_interaction_network.
                                              get_measurements(
                                                  network, time, modification,
                                                  site_combination,
                                                  replicate_combination)))))

                    elif visual_property_sub_element.get(
                            "name") == f"NODE_CUSTOMGRAPHICS_POSITION_{i + 1}":
                        visual_property_sub_element.set(
                            "default",
                            f"{('W', 'E')[i]},{('E', 'W')[i]},c,0.00,0.00",
                        )

    ET.indent(styles)
    return styles


def export(styles: ET.ElementTree, basename: str) -> None:
    """
    Exports the Cytoscape styles.

    Args:
        styles: The Cytoscape styles.
        basename: The base file name.
    """
    styles.write(
        f"{basename}.xml",
        encoding="utf-8",
        xml_declaration=True,
    )
