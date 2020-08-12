SHAPE_MAP = [{
    "tag": "discreteMappingEntry",
    "attrib": {
        "attributeValue": "P",
        "value": "ROUND_RECTANGLE"
    },
    "children": []
}, {
    "tag": "discreteMappingEntry",
    "attrib": {
        "attributeValue": "B",
        "value": "ELLIPSE"
    },
    "children": []
}, {
    "tag": "discreteMappingEntry",
    "attrib": {
        "attributeValue": "U",
        "value": "TRIANGLE"
    },
    "children": []
}, {
    "tag": "discreteMappingEntry",
    "attrib": {
        "attributeValue": "N",
        "value": "DIAMOND"
    },
    "children": []
}]

COLOR_MAP = [{
    "tag": "discreteMappingEntry",
    "attrib": {
        "attributeValue": "UP",
        "value": "#FF0000"
    },
    "children": []
}, {
    "tag": "discreteMappingEntry",
    "attrib": {
        "attributeValue": "MID_UP",
        "value": "#FF6666"
    },
    "children": []
}, {
    "tag": "discreteMappingEntry",
    "attrib": {
        "attributeValue": "MID",
        "value": "#808080"
    },
    "children": []
}, {
    "tag": "discreteMappingEntry",
    "attrib": {
        "attributeValue": "MID_DOWN",
        "value": "#6666FF"
    },
    "children": []
}, {
    "tag": "discreteMappingEntry",
    "attrib": {
        "attributeValue": "DOWN",
        "value": "#0000FF"
    },
    "children": []
}, {
    "tag": "discreteMappingEntry",
    "attrib": {
        "attributeValue": "P_UP_U_DOWN",
        "value": "#00FF00"
    },
    "children": []
}, {
    "tag": "discreteMappingEntry",
    "attrib": {
        "attributeValue": "P_DOWN_U_UP",
        "value": "#FFFF00"
    },
    "children": []
}]

PTM_TO_HEX = {
    mapping["attrib"]["attributeValue"]: mapping["attrib"]["value"]
    for mapping in COLOR_MAP
}

CHART_POSITION = {"L": "W,E,c,0.00,0.00", "R": "E,W,c,0.00,0.00"}
STYLE = {
    'tag':
        'vizmap',
    'attrib': {
        'id': "PTM_PPI_NETWORK",
        'documentVersion': '3.0'
    },
    'children': [{
        'tag':
            'visualStyle',
        'attrib': {
            'name': 'default'
        },
        'children': [{
            'tag':
                'network',
            'attrib': {},
            'children': [{
                'tag': 'visualProperty',
                'attrib': {
                    'default': 'true',
                    'name': 'NETWORK_EDGE_SELECTION'
                },
                'children': []
            }, {
                'tag': 'visualProperty',
                'attrib': {
                    'default': '0.0',
                    'name': 'NETWORK_CENTER_Y_LOCATION'
                },
                'children': []
            }, {
                'tag': 'visualProperty',
                'attrib': {
                    'default': '0.0',
                    'name': 'NETWORK_CENTER_X_LOCATION'
                },
                'children': []
            }, {
                'tag': 'visualProperty',
                'attrib': {
                    'default': '1.0',
                    'name': 'NETWORK_SCALE_FACTOR'
                },
                'children': []
            }, {
                'tag': 'visualProperty',
                'attrib': {
                    'default': '',
                    'name': 'NETWORK_TITLE'
                },
                'children': []
            }, {
                'tag': 'visualProperty',
                'attrib': {
                    'default': '0.0',
                    'name': 'NETWORK_DEPTH'
                },
                'children': []
            }, {
                'tag': 'visualProperty',
                'attrib': {
                    'default': '#FFFFFF',
                    'name': 'NETWORK_BACKGROUND_PAINT'
                },
                'children': []
            }, {
                'tag': 'visualProperty',
                'attrib': {
                    'default': '550.0',
                    'name': 'NETWORK_WIDTH'
                },
                'children': []
            }, {
                'tag': 'visualProperty',
                'attrib': {
                    'default': '550.0',
                    'name': 'NETWORK_SIZE'
                },
                'children': []
            }, {
                'tag': 'visualProperty',
                'attrib': {
                    'default': '0.0',
                    'name': 'NETWORK_CENTER_Z_LOCATION'
                },
                'children': []
            }, {
                'tag': 'visualProperty',
                'attrib': {
                    'default': 'true',
                    'name': 'NETWORK_NODE_SELECTION'
                },
                'children': []
            }, {
                'tag': 'visualProperty',
                'attrib': {
                    'default': '400.0',
                    'name': 'NETWORK_HEIGHT'
                },
                'children': []
            }]
        }, {
            'tag':
                'node',
            'attrib': {},
            'children': [{
                'tag': 'dependency',
                'attrib': {
                    'value': 'true',
                    'name': 'nodeCustomGraphicsSizeSync'
                },
                'children': []
            }, {
                'tag': 'dependency',
                'attrib': {
                    'value': 'false',
                    'name': 'nodeSizeLocked'
                },
                'children': []
            }, {
                'tag': 'visualProperty',
                'attrib': {
                    'default': 'C,C,c,0.00,0.00',
                    'name': 'NODE_CUSTOMGRAPHICS_POSITION_5'
                },
                'children': []
            }, {
                'tag': 'visualProperty',
                'attrib': {
                    'default': 'C,C,c,0.00,0.00',
                    'name': 'NODE_CUSTOMGRAPHICS_POSITION_1'
                },
                'children': []
            }, {
                'tag': 'visualProperty',
                'attrib': {
                    'default':
                        'org.cytoscape.ding.customgraphics.NullCustomGraphics,0,[ Remove Graphics ],',
                    'name':
                        'NODE_CUSTOMGRAPHICS_2'
                },
                'children': []
            }, {
                'tag': 'visualProperty',
                'attrib': {
                    'default':
                        'DefaultVisualizableVisualProperty(id=NODE_CUSTOMPAINT_6, name=Node Custom Paint 6)',
                    'name':
                        'NODE_CUSTOMPAINT_6'
                },
                'children': []
            }, {
                'tag': 'visualProperty',
                'attrib': {
                    'default': 'true',
                    'name': 'NODE_VISIBLE'
                },
                'children': []
            }, {
                'tag': 'visualProperty',
                'attrib': {
                    'default': '50.0',
                    'name': 'NODE_CUSTOMGRAPHICS_SIZE_9'
                },
                'children': []
            }, {
                'tag': 'visualProperty',
                'attrib': {
                    'default': '35.0',
                    'name': 'NODE_SIZE'
                },
                'children': []
            }, {
                'tag': 'visualProperty',
                'attrib': {
                    'default': '50.0',
                    'name': 'NODE_CUSTOMGRAPHICS_SIZE_7'
                },
                'children': []
            }, {
                'tag': 'visualProperty',
                'attrib': {
                    'default': '10.0',
                    'name': 'COMPOUND_NODE_PADDING'
                },
                'children': []
            }, {
                'tag': 'visualProperty',
                'attrib': {
                    'default': 'ROUND_RECTANGLE',
                    'name': 'COMPOUND_NODE_SHAPE'
                },
                'children': []
            }, {
                'tag': 'visualProperty',
                'attrib': {
                    'default':
                        'org.cytoscape.ding.customgraphics.NullCustomGraphics,0,[ Remove Graphics ],',
                    'name':
                        'NODE_CUSTOMGRAPHICS_5'
                },
                'children': []
            }, {
                'tag': 'visualProperty',
                'attrib': {
                    'default': '#CCCCCC',
                    'name': 'NODE_BORDER_PAINT'
                },
                'children': []
            }, {
                'tag': 'visualProperty',
                'attrib': {
                    'default': 'C,C,c,0.00,0.00',
                    'name': 'NODE_CUSTOMGRAPHICS_POSITION_8'
                },
                'children': []
            }, {
                'tag': 'visualProperty',
                'attrib': {
                    'default': '#1E90FF',
                    'name': 'NODE_PAINT'
                },
                'children': []
            }, {
                'tag': 'visualProperty',
                'attrib': {
                    'default': 'false',
                    'name': 'NODE_SELECTED'
                },
                'children': []
            }, {
                'tag': 'visualProperty',
                'attrib': {
                    'default': '0.0',
                    'name': 'NODE_Y_LOCATION'
                },
                'children': []
            }, {
                'tag': 'visualProperty',
                'attrib': {
                    'default': '255',
                    'name': 'NODE_TRANSPARENCY'
                },
                'children': []
            }, {
                'tag': 'visualProperty',
                'attrib': {
                    'default':
                        'DefaultVisualizableVisualProperty(id=NODE_CUSTOMPAINT_5, name=Node Custom Paint 5)',
                    'name':
                        'NODE_CUSTOMPAINT_5'
                },
                'children': []
            }, {
                'tag': 'visualProperty',
                'attrib': {
                    'default': '',
                    'name': 'NODE_TOOLTIP'
                },
                'children': []
            }, {
                'tag': 'visualProperty',
                'attrib': {
                    'default':
                        'org.cytoscape.ding.customgraphics.NullCustomGraphics,0,[ Remove Graphics ],',
                    'name':
                        'NODE_CUSTOMGRAPHICS_7'
                },
                'children': []
            }, {
                'tag': 'visualProperty',
                'attrib': {
                    'default': '12',
                    'name': 'NODE_LABEL_FONT_SIZE'
                },
                'children': []
            }, {
                'tag': 'visualProperty',
                'attrib': {
                    'default': 'C,C,c,0.00,0.00',
                    'name': 'NODE_CUSTOMGRAPHICS_POSITION_2'
                },
                'children': []
            }, {
                'tag': 'visualProperty',
                'attrib': {
                    'default':
                        'org.cytoscape.ding.customgraphics.NullCustomGraphics,0,[ Remove Graphics ],',
                    'name':
                        'NODE_CUSTOMGRAPHICS_4'
                },
                'children': []
            }, {
                'tag': 'visualProperty',
                'attrib': {
                    'default': 'true',
                    'name': 'NODE_NESTED_NETWORK_IMAGE_VISIBLE'
                },
                'children': []
            }, {
                'tag': 'visualProperty',
                'attrib': {
                    'default': '255',
                    'name': 'NODE_LABEL_TRANSPARENCY'
                },
                'children': []
            }, {
                'tag': 'visualProperty',
                'attrib': {
                    'default':
                        'DefaultVisualizableVisualProperty(id=NODE_CUSTOMPAINT_3, name=Node Custom Paint 3)',
                    'name':
                        'NODE_CUSTOMPAINT_3'
                },
                'children': []
            }, {
                'tag': 'visualProperty',
                'attrib': {
                    'default':
                        'DefaultVisualizableVisualProperty(id=NODE_CUSTOMPAINT_2, name=Node Custom Paint 2)',
                    'name':
                        'NODE_CUSTOMPAINT_2'
                },
                'children': []
            }, {
                'tag': 'visualProperty',
                'attrib': {
                    'default': '75.0',
                    'name': 'NODE_WIDTH'
                },
                'children': []
            }, {
                'tag': 'visualProperty',
                'attrib': {
                    'default': '0.0',
                    'name': 'NODE_DEPTH'
                },
                'children': []
            }, {
                'tag': 'visualProperty',
                'attrib': {
                    'default': '50.0',
                    'name': 'NODE_CUSTOMGRAPHICS_SIZE_5'
                },
                'children': []
            }, {
                'tag': 'visualProperty',
                'attrib': {
                    'default': 'C,C,c,0.00,0.00',
                    'name': 'NODE_CUSTOMGRAPHICS_POSITION_7'
                },
                'children': []
            }, {
                'tag': 'visualProperty',
                'attrib': {
                    'default': '0.0',
                    'name': 'NODE_X_LOCATION'
                },
                'children': []
            }, {
                'tag': 'visualProperty',
                'attrib': {
                    'default':
                        'org.cytoscape.ding.customgraphics.NullCustomGraphics,0,[ Remove Graphics ],',
                    'name':
                        'NODE_CUSTOMGRAPHICS_3'
                },
                'children': []
            }, {
                'tag': 'visualProperty',
                'attrib': {
                    'default': 'C,C,c,0.00,0.00',
                    'name': 'NODE_LABEL_POSITION'
                },
                'children': []
            }, {
                'tag': 'visualProperty',
                'attrib': {
                    'default':
                        'org.cytoscape.ding.customgraphics.NullCustomGraphics,0,[ Remove Graphics ],',
                    'name':
                        'NODE_CUSTOMGRAPHICS_1'
                },
                'children': []
            }, {
                'tag': 'visualProperty',
                'attrib': {
                    'default': '#FFFF00',
                    'name': 'NODE_SELECTED_PAINT'
                },
                'children': []
            }, {
                'tag': 'visualProperty',
                'attrib': {
                    'default':
                        'DefaultVisualizableVisualProperty(id=NODE_CUSTOMPAINT_8, name=Node Custom Paint 8)',
                    'name':
                        'NODE_CUSTOMPAINT_8'
                },
                'children': []
            }, {
                'tag': 'visualProperty',
                'attrib': {
                    'default': 'SansSerif.plain,plain,12',
                    'name': 'NODE_LABEL_FONT_FACE'
                },
                'children': []
            }, {
                'tag': 'visualProperty',
                'attrib': {
                    'default': '35.0',
                    'name': 'NODE_HEIGHT'
                },
                'children': []
            }, {
                'tag': 'visualProperty',
                'attrib': {
                    'default': '200.0',
                    'name': 'NODE_LABEL_WIDTH'
                },
                'children': []
            }, {
                'tag': 'visualProperty',
                'attrib': {
                    'default': '#000000',
                    'name': 'NODE_LABEL_COLOR'
                },
                'children': []
            }, {
                'tag': 'visualProperty',
                'attrib': {
                    'default': '0.0',
                    'name': 'NODE_BORDER_WIDTH'
                },
                'children': []
            }, {
                'tag': 'visualProperty',
                'attrib': {
                    'default':
                        'DefaultVisualizableVisualProperty(id=NODE_CUSTOMPAINT_9, name=Node Custom Paint 9)',
                    'name':
                        'NODE_CUSTOMPAINT_9'
                },
                'children': []
            }, {
                'tag': 'visualProperty',
                'attrib': {
                    'default':
                        'org.cytoscape.ding.customgraphics.NullCustomGraphics,0,[ Remove Graphics ],',
                    'name':
                        'NODE_CUSTOMGRAPHICS_9'
                },
                'children': []
            }, {
                'tag': 'visualProperty',
                'attrib': {
                    'default': '50.0',
                    'name': 'NODE_CUSTOMGRAPHICS_SIZE_4'
                },
                'children': []
            }, {
                'tag': 'visualProperty',
                'attrib': {
                    'default':
                        'DefaultVisualizableVisualProperty(id=NODE_CUSTOMPAINT_7, name=Node Custom Paint 7)',
                    'name':
                        'NODE_CUSTOMPAINT_7'
                },
                'children': []
            }, {
                'tag': 'visualProperty',
                'attrib': {
                    'default': '50.0',
                    'name': 'NODE_CUSTOMGRAPHICS_SIZE_8'
                },
                'children': []
            }, {
                'tag': 'visualProperty',
                'attrib': {
                    'default':
                        'DefaultVisualizableVisualProperty(id=NODE_CUSTOMPAINT_1, name=Node Custom Paint 1)',
                    'name':
                        'NODE_CUSTOMPAINT_1'
                },
                'children': []
            }, {
                'tag': 'visualProperty',
                'attrib': {
                    'default': '50.0',
                    'name': 'NODE_CUSTOMGRAPHICS_SIZE_3'
                },
                'children': []
            }, {
                'tag': 'visualProperty',
                'attrib': {
                    'default': '50.0',
                    'name': 'NODE_CUSTOMGRAPHICS_SIZE_1'
                },
                'children': []
            }, {
                'tag': 'visualProperty',
                'attrib': {
                    'default': '50.0',
                    'name': 'NODE_CUSTOMGRAPHICS_SIZE_2'
                },
                'children': []
            }, {
                'tag': 'visualProperty',
                'attrib': {
                    'default': 'C,C,c,0.00,0.00',
                    'name': 'NODE_CUSTOMGRAPHICS_POSITION_9'
                },
                'children': []
            }, {
                'tag': 'visualProperty',
                'attrib': {
                    'default':
                        'org.cytoscape.ding.customgraphics.NullCustomGraphics,0,[ Remove Graphics ],',
                    'name':
                        'NODE_CUSTOMGRAPHICS_6'
                },
                'children': []
            }, {
                'tag': 'visualProperty',
                'attrib': {
                    'default': '0.0',
                    'name': 'NODE_Z_LOCATION'
                },
                'children': []
            }, {
                'tag': 'visualProperty',
                'attrib': {
                    'default': 'SOLID',
                    'name': 'NODE_BORDER_STROKE'
                },
                'children': []
            }, {
                'tag': 'visualProperty',
                'attrib': {
                    'default': 'C,C,c,0.00,0.00',
                    'name': 'NODE_CUSTOMGRAPHICS_POSITION_3'
                },
                'children': []
            }, {
                'tag': 'visualProperty',
                'attrib': {
                    'default': 'C,C,c,0.00,0.00',
                    'name': 'NODE_CUSTOMGRAPHICS_POSITION_4'
                },
                'children': []
            }, {
                'tag': 'visualProperty',
                'attrib': {
                    'default': '50.0',
                    'name': 'NODE_CUSTOMGRAPHICS_SIZE_6'
                },
                'children': []
            }, {
                'tag':
                    'visualProperty',
                'attrib': {
                    'default': '',
                    'name': 'NODE_LABEL'
                },
                'children': [{
                    'tag': 'passthroughMapping',
                    'attrib': {
                        'attributeName': 'name',
                        'attributeType': 'string'
                    },
                    'children': []
                }]
            }, {
                'tag': 'visualProperty',
                'attrib': {
                    'default': 'C,C,c,0.00,0.00',
                    'name': 'NODE_CUSTOMGRAPHICS_POSITION_6'
                },
                'children': []
            }, {
                'tag': 'visualProperty',
                'attrib': {
                    'default': '255',
                    'name': 'NODE_BORDER_TRANSPARENCY'
                },
                'children': []
            }, {
                'tag': 'visualProperty',
                'attrib': {
                    'default': '#89D0F5',
                    'name': 'NODE_FILL_COLOR'
                },
                'children': []
            }, {
                'tag': 'visualProperty',
                'attrib': {
                    'default': 'ROUND_RECTANGLE',
                    'name': 'NODE_SHAPE'
                },
                'children': []
            }, {
                'tag': 'visualProperty',
                'attrib': {
                    'default':
                        'DefaultVisualizableVisualProperty(id=NODE_CUSTOMPAINT_4, name=Node Custom Paint 4)',
                    'name':
                        'NODE_CUSTOMPAINT_4'
                },
                'children': []
            }, {
                'tag': 'visualProperty',
                'attrib': {
                    'default':
                        'org.cytoscape.ding.customgraphics.NullCustomGraphics,0,[ Remove Graphics ],',
                    'name':
                        'NODE_CUSTOMGRAPHICS_8'
                },
                'children': []
            }]
        }, {
            'tag':
                'edge',
            'attrib': {},
            'children': [{
                'tag': 'dependency',
                'attrib': {
                    'value': 'false',
                    'name': 'arrowColorMatchesEdge'
                },
                'children': []
            }, {
                'tag': 'visualProperty',
                'attrib': {
                    'default': '255',
                    'name': 'EDGE_TRANSPARENCY'
                },
                'children': []
            }, {
                'tag': 'visualProperty',
                'attrib': {
                    'default': '255',
                    'name': 'EDGE_LABEL_TRANSPARENCY'
                },
                'children': []
            }, {
                'tag': 'visualProperty',
                'attrib': {
                    'default': '',
                    'name': 'EDGE_TOOLTIP'
                },
                'children': []
            }, {
                'tag': 'visualProperty',
                'attrib': {
                    'default': '#323232',
                    'name': 'EDGE_PAINT'
                },
                'children': []
            }, {
                'tag': 'visualProperty',
                'attrib': {
                    'default': '#000000',
                    'name': 'EDGE_TARGET_ARROW_UNSELECTED_PAINT'
                },
                'children': []
            }, {
                'tag': 'visualProperty',
                'attrib': {
                    'default': '',
                    'name': 'EDGE_LABEL'
                },
                'children': []
            }, {
                'tag': 'visualProperty',
                'attrib': {
                    'default': 'false',
                    'name': 'EDGE_SELECTED'
                },
                'children': []
            }, {
                'tag': 'visualProperty',
                'attrib': {
                    'default': '#FF0000',
                    'name': 'EDGE_SELECTED_PAINT'
                },
                'children': []
            }, {
                'tag': 'visualProperty',
                'attrib': {
                    'default': '#FFFF00',
                    'name': 'EDGE_SOURCE_ARROW_SELECTED_PAINT'
                },
                'children': []
            }, {
                'tag': 'visualProperty',
                'attrib': {
                    'default': '#FFFF00',
                    'name': 'EDGE_TARGET_ARROW_SELECTED_PAINT'
                },
                'children': []
            }, {
                'tag': 'visualProperty',
                'attrib': {
                    'default': '2.0',
                    'name': 'EDGE_WIDTH'
                },
                'children': []
            }, {
                'tag': 'visualProperty',
                'attrib': {
                    'default': 'SOLID',
                    'name': 'EDGE_LINE_TYPE'
                },
                'children': []
            }, {
                'tag': 'visualProperty',
                'attrib': {
                    'default': '200.0',
                    'name': 'EDGE_LABEL_WIDTH'
                },
                'children': []
            }, {
                'tag': 'visualProperty',
                'attrib': {
                    'default': '#000000',
                    'name': 'EDGE_LABEL_COLOR'
                },
                'children': []
            }, {
                'tag': 'visualProperty',
                'attrib': {
                    'default': '10',
                    'name': 'EDGE_LABEL_FONT_SIZE'
                },
                'children': []
            }, {
                'tag': 'visualProperty',
                'attrib': {
                    'default': '',
                    'name': 'EDGE_BEND'
                },
                'children': []
            }, {
                'tag': 'visualProperty',
                'attrib': {
                    'default': 'NONE',
                    'name': 'EDGE_TARGET_ARROW_SHAPE'
                },
                'children': []
            }, {
                'tag': 'visualProperty',
                'attrib': {
                    'default': '#000000',
                    'name': 'EDGE_SOURCE_ARROW_UNSELECTED_PAINT'
                },
                'children': []
            }, {
                'tag': 'visualProperty',
                'attrib': {
                    'default': 'Dialog.plain,plain,10',
                    'name': 'EDGE_LABEL_FONT_FACE'
                },
                'children': []
            }, {
                'tag': 'visualProperty',
                'attrib': {
                    'default': '#848484',
                    'name': 'EDGE_STROKE_UNSELECTED_PAINT'
                },
                'children': []
            }, {
                'tag': 'visualProperty',
                'attrib': {
                    'default': '#FF0000',
                    'name': 'EDGE_STROKE_SELECTED_PAINT'
                },
                'children': []
            }, {
                'tag': 'visualProperty',
                'attrib': {
                    'default': 'true',
                    'name': 'EDGE_CURVED'
                },
                'children': []
            }, {
                'tag': 'visualProperty',
                'attrib': {
                    'default': 'NONE',
                    'name': 'EDGE_SOURCE_ARROW_SHAPE'
                },
                'children': []
            }, {
                'tag': 'visualProperty',
                'attrib': {
                    'default': 'true',
                    'name': 'EDGE_VISIBLE'
                },
                'children': []
            }, {
                'tag': 'visualProperty',
                'attrib': {
                    'default': '#404040',
                    'name': 'EDGE_UNSELECTED_PAINT'
                },
                'children': []
            }]
        }]
    }]
}
