"""Cytoscape style configurations for a protein-protein interaction network."""
COMPONENTS = [
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
                    "default":
                        "org.cytoscape.ding.customgraphics.NullCustomGraphics,0,[ Remove Graphics ],"
                },
                "NODE_CUSTOMGRAPHICS_2": {
                    "default":
                        "org.cytoscape.ding.customgraphics.NullCustomGraphics,0,[ Remove Graphics ],"
                },
                "NODE_CUSTOMGRAPHICS_3": {
                    "default":
                        "org.cytoscape.ding.customgraphics.NullCustomGraphics,0,[ Remove Graphics ],"
                },
                "NODE_CUSTOMGRAPHICS_4": {
                    "default":
                        "org.cytoscape.ding.customgraphics.NullCustomGraphics,0,[ Remove Graphics ],"
                },
                "NODE_CUSTOMGRAPHICS_5": {
                    "default":
                        "org.cytoscape.ding.customgraphics.NullCustomGraphics,0,[ Remove Graphics ],"
                },
                "NODE_CUSTOMGRAPHICS_6": {
                    "default":
                        "org.cytoscape.ding.customgraphics.NullCustomGraphics,0,[ Remove Graphics ],"
                },
                "NODE_CUSTOMGRAPHICS_7": {
                    "default":
                        "org.cytoscape.ding.customgraphics.NullCustomGraphics,0,[ Remove Graphics ],"
                },
                "NODE_CUSTOMGRAPHICS_8": {
                    "default":
                        "org.cytoscape.ding.customgraphics.NullCustomGraphics,0,[ Remove Graphics ],"
                },
                "NODE_CUSTOMGRAPHICS_9": {
                    "default":
                        "org.cytoscape.ding.customgraphics.NullCustomGraphics,0,[ Remove Graphics ],"
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
                        "DefaultVisualizableVisualProperty(id=NODE_CUSTOMPAINT_1, name=Node Custom Paint 1)"
                },
                "NODE_CUSTOMPAINT_2": {
                    "default":
                        "DefaultVisualizableVisualProperty(id=NODE_CUSTOMPAINT_2, name=Node Custom Paint 2)"
                },
                "NODE_CUSTOMPAINT_3": {
                    "default":
                        "DefaultVisualizableVisualProperty(id=NODE_CUSTOMPAINT_3, name=Node Custom Paint 3)"
                },
                "NODE_CUSTOMPAINT_4": {
                    "default":
                        "DefaultVisualizableVisualProperty(id=NODE_CUSTOMPAINT_4, name=Node Custom Paint 4)"
                },
                "NODE_CUSTOMPAINT_5": {
                    "default":
                        "DefaultVisualizableVisualProperty(id=NODE_CUSTOMPAINT_5, name=Node Custom Paint 5)"
                },
                "NODE_CUSTOMPAINT_6": {
                    "default":
                        "DefaultVisualizableVisualProperty(id=NODE_CUSTOMPAINT_6, name=Node Custom Paint 6)"
                },
                "NODE_CUSTOMPAINT_7": {
                    "default":
                        "DefaultVisualizableVisualProperty(id=NODE_CUSTOMPAINT_7, name=Node Custom Paint 7)"
                },
                "NODE_CUSTOMPAINT_8": {
                    "default":
                        "DefaultVisualizableVisualProperty(id=NODE_CUSTOMPAINT_8, name=Node Custom Paint 8)"
                },
                "NODE_CUSTOMPAINT_9": {
                    "default":
                        "DefaultVisualizableVisualProperty(id=NODE_CUSTOMPAINT_9, name=Node Custom Paint 9)"
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
                        "attributeName": "confidence",
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
                    "default":
                        "org.cytoscape.ding.customgraphics.NullCustomGraphics,0,[ Remove Graphics ],"
                },
                "NODE_CUSTOMGRAPHICS_2": {
                    "default":
                        "org.cytoscape.ding.customgraphics.NullCustomGraphics,0,[ Remove Graphics ],"
                },
                "NODE_CUSTOMGRAPHICS_3": {
                    "default":
                        "org.cytoscape.ding.customgraphics.NullCustomGraphics,0,[ Remove Graphics ],"
                },
                "NODE_CUSTOMGRAPHICS_4": {
                    "default":
                        "org.cytoscape.ding.customgraphics.NullCustomGraphics,0,[ Remove Graphics ],"
                },
                "NODE_CUSTOMGRAPHICS_5": {
                    "default":
                        "org.cytoscape.ding.customgraphics.NullCustomGraphics,0,[ Remove Graphics ],"
                },
                "NODE_CUSTOMGRAPHICS_6": {
                    "default":
                        "org.cytoscape.ding.customgraphics.NullCustomGraphics,0,[ Remove Graphics ],"
                },
                "NODE_CUSTOMGRAPHICS_7": {
                    "default":
                        "org.cytoscape.ding.customgraphics.NullCustomGraphics,0,[ Remove Graphics ],"
                },
                "NODE_CUSTOMGRAPHICS_8": {
                    "default":
                        "org.cytoscape.ding.customgraphics.NullCustomGraphics,0,[ Remove Graphics ],"
                },
                "NODE_CUSTOMGRAPHICS_9": {
                    "default":
                        "org.cytoscape.ding.customgraphics.NullCustomGraphics,0,[ Remove Graphics ],"
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
                        "DefaultVisualizableVisualProperty(id=NODE_CUSTOMPAINT_1, name=Node Custom Paint 1)"
                },
                "NODE_CUSTOMPAINT_2": {
                    "default":
                        "DefaultVisualizableVisualProperty(id=NODE_CUSTOMPAINT_2, name=Node Custom Paint 2)"
                },
                "NODE_CUSTOMPAINT_3": {
                    "default":
                        "DefaultVisualizableVisualProperty(id=NODE_CUSTOMPAINT_3, name=Node Custom Paint 3)"
                },
                "NODE_CUSTOMPAINT_4": {
                    "default":
                        "DefaultVisualizableVisualProperty(id=NODE_CUSTOMPAINT_4, name=Node Custom Paint 4)"
                },
                "NODE_CUSTOMPAINT_5": {
                    "default":
                        "DefaultVisualizableVisualProperty(id=NODE_CUSTOMPAINT_5, name=Node Custom Paint 5)"
                },
                "NODE_CUSTOMPAINT_6": {
                    "default":
                        "DefaultVisualizableVisualProperty(id=NODE_CUSTOMPAINT_6, name=Node Custom Paint 6)"
                },
                "NODE_CUSTOMPAINT_7": {
                    "default":
                        "DefaultVisualizableVisualProperty(id=NODE_CUSTOMPAINT_7, name=Node Custom Paint 7)"
                },
                "NODE_CUSTOMPAINT_8": {
                    "default":
                        "DefaultVisualizableVisualProperty(id=NODE_CUSTOMPAINT_8, name=Node Custom Paint 8)"
                },
                "NODE_CUSTOMPAINT_9": {
                    "default":
                        "DefaultVisualizableVisualProperty(id=NODE_CUSTOMPAINT_9, name=Node Custom Paint 9)"
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
                            "{modifications[0]}_up {modifications[1]}_down":
                                f"#{0:02X}{255:02X}{0:02X}",
                            "{modifications[0]}_up {modifications[1]}_mid_down":
                                f"#{0:02X}{255:02X}{0:02X}",
                            "{modifications[0]}_mid_up {modifications[1]}_down":
                                f"#{0:02X}{255:02X}{0:02X}",
                            "{modifications[0]}_mid_up {modifications[1]}_mid_down":
                                f"#{0:02X}{255:02X}{0:02X}",
                            "{modifications[0]}_down {modifications[1]}_up":
                                f"#{255:02X}{255:02X}{0:02X}",
                            "{modifications[0]}_down {modifications[1]}_mid_up":
                                f"#{255:02X}{255:02X}{0:02X}",
                            "{modifications[0]}_mid_down {modifications[1]}_up":
                                f"#{255:02X}{255:02X}{0:02X}",
                            "{modifications[0]}_mid_down {modifications[1]}_mid_up":
                                f"#{255:02X}{255:02X}{0:02X}",
                            "{modifications[0]}_up {modifications[1]}_mid":
                                f"#{255:02X}{0:02X}{0:02X}",
                            "{modifications[0]}_mid_up {modifications[1]}_mid":
                                f"#{255:02X}{102:02X}{102:02X}",
                            "{modifications[0]}_mid_down {modifications[1]}_mid":
                                f"#{102:02X}{102:02X}{255:02X}",
                            "{modifications[0]}_down {modifications[1]}_mid":
                                f"#{0:02X}{0:02X}{255:02X}",
                            "{modifications[0]}_mid {modifications[1]}_up":
                                f"#{255:02X}{0:02X}{0:02X}",
                            "{modifications[0]}_mid {modifications[1]}_mid_up":
                                f"#{255:02X}{102:02X}{102:02X}",
                            "{modifications[0]}_mid {modifications[1]}_mid_down":
                                f"#{102:02X}{102:02X}{255:02X}",
                            "{modifications[0]}_mid {modifications[1]}_down":
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
                            "{modifications[0]}": "RECTANGLE",
                            "{modifications[1]}": "TRIANGLE",
                            "{modifications[0]}{modifications[1]}": "ELLIPSE",
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
