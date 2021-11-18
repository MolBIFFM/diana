import json
import xml.etree.ElementTree as ET

from cytoscape import configuration
from networks import protein_interaction_network


def get_styles(
        network,
        bar_chart_range=(-3.0, 3.0),
        get_bar_chart_range=lambda time, ptm, bar_chart_range,
    site_combination: bar_chart_range,
        site_combination=lambda sites: max(sites, key=abs),
):
    styles = ET.ElementTree(
        ET.Element("vizmap", attrib={
            "id": "VizMap",
            "documentVersion": "3.0"
        }))

    for time in protein_interaction_network.get_times(network):
        modifications = protein_interaction_network.get_post_translational_modifications(
            network, time)
        visual_style = ET.SubElement(styles.getroot(),
                                     "visualStyle",
                                     attrib={"name": str(time)})
        tags = {}
        for tag in configuration.SETTINGS[len(modifications) - 1]:
            tags[tag] = ET.SubElement(visual_style, tag)
            for name, setting in configuration.SETTINGS[
                    len(modifications) - 1][tag]["dependency"].items():
                ET.SubElement(
                    tags[tag],
                    "dependency",
                    attrib={
                        "name": name,
                        "value": setting["value"]
                    },
                )
            for name, setting in configuration.SETTINGS[
                    len(modifications) - 1][tag]["visualProperty"].items():
                visual_property = ET.SubElement(
                    tags[tag],
                    "visualProperty",
                    attrib={
                        "name": name,
                        "default": setting["default"]
                    },
                )
                if setting.get("passthroughMapping"):
                    ET.SubElement(
                        visual_property,
                        "passthroughMapping",
                        attrib={
                            "attributeName":
                            setting["passthroughMapping"]["attributeName"],
                            "attributeType":
                            setting["passthroughMapping"]["attributeType"],
                        },
                    )
                elif setting.get("discreteMapping"):
                    discrete_mapping = ET.SubElement(
                        visual_property,
                        "discreteMapping",
                        attrib={
                            "attributeName":
                            setting["discreteMapping"]["attributeName"].format(
                                time=time),
                            "attributeType":
                            setting["discreteMapping"]["attributeType"],
                        },
                    )
                    for attribute_value, value in setting["discreteMapping"][
                            "discreteMappingEntry"].items():
                        ET.SubElement(
                            discrete_mapping,
                            "discreteMappingEntry",
                            attrib={
                                "attributeValue":
                                attribute_value.format(
                                    modifications=modifications),
                                "value":
                                value,
                            },
                        )

        for i, ptm in enumerate(modifications):
            for visual_property in visual_style.find("node").findall(
                    "visualProperty"):
                if visual_property.get(
                        "name") == "NODE_CUSTOMGRAPHICS_{}".format(i + 1):
                    visual_property.set(
                        "default",
                        get_bar_chart(
                            time,
                            ptm,
                            protein_interaction_network.get_sites(
                                network, time, ptm),
                            cy_range=get_bar_chart_range(
                                time, ptm, bar_chart_range, site_combination),
                        ),
                    )
                elif visual_property.get(
                        "name") == "NODE_CUSTOMGRAPHICS_POSITION_{}".format(i +
                                                                            1):
                    visual_property.set(
                        "default",
                        "{},{},c,0.00,0.00".format(*[("W", "E"), ("E",
                                                                  "W")][i]),
                    )

    return styles


def get_bar_chart(time, ptm, sites, cy_range=(-2.0, 2.0)):
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
        "cy_dataColumns":
        ["{} {} {}".format(time, ptm, site + 1) for site in range(sites)],
    })
    return "org.cytoscape.BarChart: {}".format(bar_chart)


def export(styles, basename, suffix=""):
    styles.write(
        "{0}{1}.xml".format(basename, suffix),
        encoding="UTF-8",
        xml_declaration=True,
    )
