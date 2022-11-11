"""
The Cytoscape style configuration for a protein-protein interaction
network.
"""

import json
import statistics
import xml.etree.ElementTree as ET
from typing import Callable, Collection, Iterable, Optional

import networkx as nx
from cytoscape import elements
from networks import protein_interaction_network

COMPONENTS: dict[str, dict[str, dict[str, dict[
    str, bool | int | str | float | dict[str, bool | int | str | float | dict[
        str, bool | int | str |
        float | dict[str, bool | int | str | float]]]]]]] = {
            "edge": {
                "dependency": {
                    "arrowColorMatchesEdge": {
                        "value": False
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
                        "default": "NONE"
                    },
                    "EDGE_TARGET_ARROW_UNSELECTED_PAINT": {
                        "default": f"#{0:02X}{0:02X}{0:02X}"
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
                        "default": 0.0
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
                        "default": f"#{128:02X}{128:02X}{128:02X}",
                    },
                    "NODE_HEIGHT": {
                        "default": 35.0
                    },
                    "NODE_LABEL": {
                        "default": ""
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
                        "default": "ROUND_RECTANGLE",
                    },
                    "NODE_SIZE": {
                        "default": 35.0,
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
                }
            }
        }


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
            f"{time} {modification} S{site}" for site in range(1, sites + 1)
        ],
    })
    return f"org.cytoscape.BarChart: {bar_chart}"


def get_styles(
    network: nx.Graph, node_shape_modifications: Iterable[str],
    node_color_modifications: Iterable[str],
    node_size_modification: Optional[str],
    bar_chart_modifications: Iterable[str],
    measurement_conversion: dict[str, Callable[[float, Collection[float]],
                                               float]],
    site_average: dict[str, Callable[[Iterable[float]], float]],
    replicate_average: dict[str, Callable[[Iterable[float]], float]],
    confidence_score_average: Callable[[dict[str, float]], float]
) -> ET.ElementTree:
    """
    Returns the Cytoscape styles for a protein-protein interaction network.

    Args:
        network: The protein-protein interaction network.
        node_shape_modifications: The identifiers of
            post-translational modifications to represent by node shape.
        node_color_modifications: The identifiers of
            post-translational modifications to represent by node color.
        node_size_modification: The identifier of a
            post-translational modification to represent by node size.
        bar_chart_modifications: The identifiers of site-specific
            post-translational modifications to represent by bar charts.
        measurement_conversion: Modification-specific functions to transform
            binary logarithms of measurements.
        site_average: Modification-specific functions to derive
            protein-specific measurements from their site-specific measurements.
        replicate_average: Modification-specific functions to derive
            site-specific measurements from their replicate measurements.
        confidence_score_average: The function to derive an edge-specific
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

    max_edge_score = confidence_score_average({
        database: 1.0
        for database in protein_interaction_network.get_databases(network)
    })

    for time in protein_interaction_network.get_times(network):
        visual_style_sub_element = ET.SubElement(styles.getroot(),
                                                 "visualStyle",
                                                 attrib={"name": str(time)})
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
                if isinstance(visual_property["default"],
                              (bool, int, str, float)):
                    visual_properties[component][
                        name] = elements.add_visual_property(
                            component_sub_element, name,
                            visual_property["default"])

        elements.add_continuous_mapping(
            visual_properties["edge"]["EDGE_TRANSPARENCY"], "score", "float", {
                0.0: (0, 0, 0),
                max_edge_score: (255, 255, 255)
            })

        elements.add_passthrough_mapping(
            visual_properties["node"]["NODE_LABEL"], "name", "string")
        elements.add_passthrough_mapping(
            visual_properties["node"]["NODE_TOOLTIP"], "protein", "string")

        if (node_size_modification is not None and
                protein_interaction_network.is_modification(
                    network, time, node_size_modification)):
            measurements = protein_interaction_network.get_measurements(
                network, time, node_size_modification,
                site_average.get(node_size_modification,
                                 lambda sites: max(sites, key=abs)),
                replicate_average.get(node_size_modification, statistics.mean))
            if isinstance(
                    COMPONENTS["node"]["visualProperty"]["NODE_SIZE"]
                ["default"], float):
                elements.add_continuous_mapping(
                    visual_properties["node"]["NODE_SIZE"],
                    f"{time} {node_size_modification}", "float", {
                        min(measurements):
                            (COMPONENTS["node"]["visualProperty"]["NODE_SIZE"]
                             ["default"] * min(measurements),) * 3,
                        max(measurements):
                            (COMPONENTS["node"]["visualProperty"]["NODE_SIZE"]
                             ["default"] * max(measurements),) * 3
                    })

        node_shape_modifications = [
            modification for modification in node_shape_modifications
            if protein_interaction_network.is_modification(
                network, time, modification)
        ]

        if 0 < len(node_shape_modifications) < 3:
            node_shape = {}
            for value in set(
                    network.nodes[protein][str(time)] for protein in network):
                if not isinstance(value, str) or not value:
                    continue

                if len(node_shape_modifications) == 1:
                    if node_shape_modifications[0] in value.split(" ")[::2]:
                        node_shape[value] = "RECTANGLE"
                    else:
                        if isinstance(
                                COMPONENTS["node"]["visualProperty"]
                            ["NODE_SHAPE"]["default"], str):
                            node_shape[value] = COMPONENTS["node"][
                                "visualProperty"]["NODE_SHAPE"]["default"]
                else:
                    if node_shape_modifications[0] in value.split(
                            " ")[::2] and node_shape_modifications[
                                1] in value.split(" ")[::2]:
                        node_shape[value] = "ELLIPSE"
                    else:
                        if node_shape_modifications[0] in value.split(" ")[::2]:
                            node_shape[value] = "RECTANGLE"
                        elif node_shape_modifications[1] in value.split(
                                " ")[::2]:
                            node_shape[value] = "TRIANGLE"
                        else:
                            if isinstance(
                                    COMPONENTS["node"]["visualProperty"]
                                ["NODE_SHAPE"]["default"], str):
                                node_shape[value] = COMPONENTS["node"][
                                    "visualProperty"]["NODE_SHAPE"]["default"]
        else:
            if isinstance(
                    COMPONENTS["node"]["visualProperty"]["NODE_SHAPE"]
                ["default"], str):
                node_shape = {
                    value: COMPONENTS["node"]["visualProperty"]["NODE_SHAPE"]
                    ["default"]
                    for value in set(network.nodes[protein][str(time)]
                                     for protein in network)
                }

        elements.add_discrete_mapping(visual_properties["node"]["NODE_SHAPE"],
                                      str(time), "string", node_shape)

        node_color_modifications = [
            modification for modification in node_color_modifications
            if protein_interaction_network.is_modification(
                network, time, modification)
        ]
        if 0 < len(node_color_modifications) < 3:
            node_color = {}
            for value in set(
                    network.nodes[protein][str(time)] for protein in network):
                if not isinstance(value, str) or not value:
                    continue

                if len(node_color_modifications) == 1:
                    if node_color_modifications[0] in value.split(" ")[::2]:
                        summary = value.split(" ")[value.split(" ").index(
                            node_color_modifications[0]) + 1]
                        if summary == "UP":
                            node_color[value] = f"#{255:02X}{0:02X}{0:02X}"
                        elif summary == "MID_UP":
                            node_color[value] = f"#{255:02X}{102:02X}{102:02X}"
                        elif summary == "MID_DOWN":
                            node_color[value] = f"#{102:02X}{102:02X}{255:02X}"
                        elif summary == "DOWN":
                            node_color[value] = f"#{0:02X}{0:02X}{255:02X}"
                        else:
                            node_color[value] = f"#{128:02X}{128:02X}{128:02X}"
                else:
                    if node_color_modifications[0] in value.split(
                            " ")[::2] and node_color_modifications[
                                1] in value.split(" ")[::2]:
                        summaries = [
                            value.split(" ")[value.split(" ").index(
                                node_color_modifications[i]) + 1]
                            for i in range(2)
                        ]
                        if "UP" in summaries and not ("DOWN" in summaries or
                                                      "MID_DOWN" in summaries):
                            node_color[value] = f"#{255:02X}{0:02X}{0:02X}"
                        elif "MID_UP" in summaries and not (
                                "DOWN" in summaries or "MID_DOWN" in summaries):
                            node_color[value] = f"#{255:02X}{102:02X}{102:02X}"
                        elif "MID_DOWN" in summaries and not (
                                "UP" in summaries or "MID_UP" in summaries):
                            node_color[value] = f"#{102:02X}{102:02X}{255:02X}"
                        elif "DOWN" in summaries and not (
                                "UP" in summaries or "MID_UP" in summaries):
                            node_color[value] = f"#{0:02X}{0:02X}{255:02X}"
                        elif summaries[0] in ("UP",
                                              "MID_UP") and summaries[1] in (
                                                  "DOWN", "MID_DOWN"):
                            node_color[value] = f"#{0:02X}{255:02X}{0:02X}"
                        elif summaries[0] in ("MID_DOWN",
                                              "DOWN") and summaries[1] in (
                                                  "UP", "MID_UP"):
                            node_color[value] = f"#{255:02X}{255:02X}{0:02X}"
                        else:
                            node_color[value] = f"#{128:02X}{128:02X}{128:02X}"
                    else:
                        if node_color_modifications[0] in value.split(" ")[::2]:
                            summary = value.split(" ")[value.split(" ").index(
                                node_color_modifications[0]) + 1]
                            if summary == "UP":
                                node_color[value] = f"#{255:02X}{0:02X}{0:02X}"
                            elif summary == "MID_UP":
                                node_color[
                                    value] = f"#{255:02X}{102:02X}{102:02X}"
                            elif summary == "MID_DOWN":
                                node_color[
                                    value] = f"#{102:02X}{102:02X}{255:02X}"
                            elif summary == "DOWN":
                                node_color[value] = f"#{0:02X}{0:02X}{255:02X}"
                            else:
                                node_color[
                                    value] = f"#{128:02X}{128:02X}{128:02X}"
                        elif node_color_modifications[1] in value.split(
                                " ")[::2]:
                            summary = value.split(" ")[value.split(" ").index(
                                node_color_modifications[1]) + 1]
                            if summary == "UP":
                                node_color[value] = f"#{255:02X}{0:02X}{0:02X}"
                            elif summary == "MID_UP":
                                node_color[
                                    value] = f"#{255:02X}{102:02X}{102:02X}"
                            elif summary == "MID_DOWN":
                                node_color[
                                    value] = f"#{102:02X}{102:02X}{255:02X}"
                            elif summary == "DOWN":
                                node_color[value] = f"#{0:02X}{0:02X}{255:02X}"
                            else:
                                node_color[
                                    value] = f"#{128:02X}{128:02X}{128:02X}"
                        else:
                            node_color[value] = f"#{128:02X}{128:02X}{128:02X}"
        else:
            node_color = {
                value: COMPONENTS["node"]["visualProperty"]["NODE_FILL_COLOR"]
                ["default"]
                for value in set(
                    network.nodes[protein][str(time)] for protein in network)
                if isinstance(value, str) and isinstance(
                    COMPONENTS["node"]["visualProperty"]["NODE_FILL_COLOR"]
                    ["default"], str)
            }

        elements.add_discrete_mapping(
            visual_properties["node"]["NODE_FILL_COLOR"], str(time), "string",
            node_color)

        bar_chart_modifications = [
            modification for modification in bar_chart_modifications
            if protein_interaction_network.is_modification(
                network, time, modification, proteins=False)
        ]
        for m in range(min(len(bar_chart_modifications), 2)):
            measurements = protein_interaction_network.get_measurements(
                network, time, bar_chart_modifications[m],
                site_average.get(bar_chart_modifications[m],
                                 lambda sites: max(sites, key=abs)),
                replicate_average.get(bar_chart_modifications[m],
                                      statistics.mean))

            visual_properties["node"][f"NODE_CUSTOMGRAPHICS_{m+1}"].set(
                "default",
                get_bar_chart(
                    time,
                    bar_chart_modifications[m],
                    max(
                        protein_interaction_network.get_sites(
                            network, time, bar_chart_modifications[m], protein)
                        for protein in network),
                    cy_range=(
                        measurement_conversion.get(
                            bar_chart_modifications[m],
                            lambda measurement, measurements: measurement)(
                                min(measurements),
                                protein_interaction_network.get_measurements(
                                    network, time, bar_chart_modifications[m],
                                    site_average.get(
                                        bar_chart_modifications[m],
                                        lambda sites: max(sites, key=abs)),
                                    replicate_average.get(
                                        bar_chart_modifications[m],
                                        statistics.mean))),
                        measurement_conversion.get(
                            bar_chart_modifications[m],
                            lambda measurement, measurements: measurement)(
                                max(measurements),
                                protein_interaction_network.get_measurements(
                                    network, time, bar_chart_modifications[m],
                                    site_average.get(
                                        bar_chart_modifications[m],
                                        lambda sites: max(sites, key=abs)),
                                    replicate_average.get(
                                        bar_chart_modifications[m],
                                        statistics.mean))))))

            visual_properties["node"][
                f"NODE_CUSTOMGRAPHICS_POSITION_{m + 1}"].set(
                    "default",
                    f"{('W', 'E')[m]},{('E', 'W')[m]},c,0.00,0.00",
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
