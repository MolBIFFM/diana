import json
import xml.etree.ElementTree as ET

from cytoscape.configuration import protein_protein_interaction_network as protein_protein_interaction_network_style
from cytoscape.configuration import pathway_network as pathway_network_style
from cytoscape.configuration import ontology_network as ontology_network_style
from networks import protein_protein_interaction_network, pathway_network, ontology_network


def get_bar_chart(time, modification, sites, cy_range=(-2.0, 2.0)):
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
            "{} {} {}".format(time, modification, site + 1)
            for site in range(sites)
        ],
    })
    return "org.cytoscape.BarChart: {}".format(bar_chart)


def get_protein_protein_interaction_network_styles(
    network,
    bar_chart_range=(-3.0, 3.0),
    get_change=lambda network, change, time, modification, site_combination:
    change,
    site_combination=lambda sites: max(sites, key=abs),
    confidence_score_combination=lambda confidence_scores: float(
        bool(confidence_scores))):
    styles = ET.ElementTree(
        ET.Element("vizmap", attrib={
            "id": "VizMap",
            "documentVersion": "3.0"
        }))

    for time in protein_protein_interaction_network.get_times(network):
        modifications = protein_protein_interaction_network.get_post_translational_modifications(
            network, time)
        visual_style_sub_element = ET.SubElement(styles.getroot(),
                                                 "visualStyle",
                                                 attrib={"name": str(time)})

        for component in protein_protein_interaction_network_style.COMPONENTS[
                len(modifications) - 1]:
            component_sub_element = ET.SubElement(visual_style_sub_element,
                                                  component)

            for name, dependency in protein_protein_interaction_network_style.COMPONENTS[
                    len(modifications) - 1][component]["dependency"].items():
                ET.SubElement(
                    component_sub_element,
                    "dependency",
                    attrib={
                        "name": name,
                        "value": dependency["value"]
                    },
                )

            for name, visual_property in protein_protein_interaction_network_style.COMPONENTS[
                    len(modifications) -
                    1][component]["visualProperty"].items():
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
                                "attrValue": key,
                                "equalValue": values["equalValue"],
                                "greaterValue": values["greaterValue"],
                                "lesserValue": values["lesserValue"]
                            },
                        )

                elif visual_property.get("discreteMapping"):
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

                    for key, value in visual_property["discreteMapping"][
                            "discreteMappingEntry"].items():
                        ET.SubElement(
                            discrete_mapping_sub_element,
                            "discreteMappingEntry",
                            attrib={
                                "attributeValue":
                                key.format(modifications=modifications),
                                "value":
                                value,
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

        visual_property_sub_elements = visual_style_sub_element.find(
            "node").findall("visualProperty")

        for i, modification in enumerate(modifications):
            if i < 2:
                for visual_property_sub_element in visual_property_sub_elements:
                    if visual_property_sub_element.get(
                            "name") == "NODE_CUSTOMGRAPHICS_{}".format(i + 1):
                        visual_property_sub_element.set(
                            "default",
                            get_bar_chart(
                                time,
                                modification,
                                protein_protein_interaction_network.get_sites(
                                    network, time, modification),
                                cy_range=(
                                    get_change(
                                        protein_protein_interaction_network,
                                        bar_chart_range[0], time, modification,
                                        site_combination),
                                    get_change(
                                        protein_protein_interaction_network,
                                        bar_chart_range[1], time, modification,
                                        site_combination))))

                    elif visual_property_sub_element.get(
                            "name"
                    ) == "NODE_CUSTOMGRAPHICS_POSITION_{}".format(i + 1):
                        visual_property_sub_element.set(
                            "default",
                            "{},{},c,0.00,0.00".format(*[("W",
                                                          "E"), ("E",
                                                                 "W")][i]),
                        )

    styles.getroot().tail = "\n"
    ET.indent(styles)

    return styles


def get_pathway_network_style(network):
    style = ET.ElementTree(
        ET.Element("vizmap", attrib={
            "id": "VizMap",
            "documentVersion": "3.0"
        }))

    visual_style_sub_element = ET.SubElement(style.getroot(),
                                             "visualStyle",
                                             attrib={"name": "pathway"})

    for component in pathway_network_style.COMPONENTS:
        component_sub_element = ET.SubElement(visual_style_sub_element,
                                              component)

        for name, dependency in pathway_network_style.COMPONENTS[component][
                "dependency"].items():
            ET.SubElement(
                component_sub_element,
                "dependency",
                attrib={
                    "name": name,
                    "value": dependency["value"]
                },
            )

        for name, visual_property in pathway_network_style.COMPONENTS[
                component]["visualProperty"].items():
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
                        visual_property["continuousMapping"]["attributeName"],
                        "attributeType":
                        visual_property["continuousMapping"]["attributeType"]
                    })

                for key, values in visual_property["continuousMapping"][
                        "continuousMappingPoint"].items():
                    ET.SubElement(
                        continuous_mapping_sub_element,
                        "continuousMappingPoint",
                        attrib={
                            "attrValue":
                            key.format(max=max(
                                pathway_network.get_pathway_sizes(
                                    network).values())),
                            "equalValue":
                            values["equalValue"],
                            "greaterValue":
                            values["greaterValue"],
                            "lesserValue":
                            values["lesserValue"]
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
                        visual_property["passthroughMapping"]["attributeName"],
                        "attributeType":
                        visual_property["passthroughMapping"]["attributeType"],
                    },
                )

    style.getroot().tail = "\n"
    ET.indent(style)

    return style


def get_ontology_network_style(network):
    style = ET.ElementTree(
        ET.Element("vizmap", attrib={
            "id": "VizMap",
            "documentVersion": "3.0"
        }))

    visual_style_sub_element = ET.SubElement(style.getroot(),
                                             "visualStyle",
                                             attrib={"name": "ontology"})

    for component in ontology_network_style.COMPONENTS:
        component_sub_element = ET.SubElement(visual_style_sub_element,
                                              component)

        for name, dependency in ontology_network_style.COMPONENTS[component][
                "dependency"].items():
            ET.SubElement(
                component_sub_element,
                "dependency",
                attrib={
                    "name": name,
                    "value": dependency["value"]
                },
            )

        for name, visual_property in ontology_network_style.COMPONENTS[
                component]["visualProperty"].items():
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
                        visual_property["continuousMapping"]["attributeName"],
                        "attributeType":
                        visual_property["continuousMapping"]["attributeType"]
                    })

                for key, values in visual_property["continuousMapping"][
                        "continuousMappingPoint"].items():
                    ET.SubElement(
                        continuous_mapping_sub_element,
                        "continuousMappingPoint",
                        attrib={
                            "attrValue":
                            key.format(max=max(
                                ontology_network.get_term_sizes(
                                    network).values())),
                            "equalValue":
                            values["equalValue"],
                            "greaterValue":
                            values["greaterValue"],
                            "lesserValue":
                            values["lesserValue"]
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
                        visual_property["passthroughMapping"]["attributeName"],
                        "attributeType":
                        visual_property["passthroughMapping"]["attributeType"],
                    },
                )

    style.getroot().tail = "\n"
    ET.indent(style)

    return style


def export(styles, basename, suffix=""):
    styles.write(
        "{0}{1}.xml".format(basename, suffix),
        encoding="utf-8",
        xml_declaration=True,
    )
