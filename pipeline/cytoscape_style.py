import xml.etree.ElementTree as ET
import json

import pipeline.configuration.cytoscape_style


class CytoscapeStyle(ET.ElementTree):
    def __init__(self, ppi_network, bar_chart_range=(-3.0, 3.0)):
        super(CytoscapeStyle, self).__init__(
            ET.Element("vizmap",
                       attrib={
                           "id": "VizMap",
                           "documentVersion": "3.0"
                       }))

        NODE_CUSTOMGRAPHICS_POSITION = [
            "{},{},c,0.00,0.00".format(*anchor_points)
            for anchor_points in (("W", "E"), ("E", "W"), ("S", "N"),
                                  ("N", "S"), ("SW", "N"), ("NE", "N"),
                                  ("SE", "N"), ("NW", "N"), ("C", "C"))
        ]

        for time in ppi_network.times():
            visual_style = ET.SubElement(self.getroot(),
                                         "visualStyle",
                                         attrib={"name": str(time)})
            tags = {}
            for tag in pipeline.configuration.cytoscape_style.SETTINGS:
                tags[tag] = ET.SubElement(visual_style, tag)
                for name, setting in pipeline.configuration.cytoscape_style.SETTINGS[
                        tag]["dependency"].items():
                    ET.SubElement(tags[tag],
                                  "dependency",
                                  attrib={
                                      "name": name,
                                      "value": setting["value"]
                                  })
                for name, setting in pipeline.configuration.cytoscape_style.SETTINGS[
                        tag]["visualProperty"].items():
                    visual_property = ET.SubElement(tags[tag],
                                                    "visualProperty",
                                                    attrib={
                                                        "name":
                                                        name,
                                                        "default":
                                                        setting["default"]
                                                    })
                    if setting.get("passthroughMapping"):
                        ET.SubElement(
                            visual_property,
                            "passthroughMapping",
                            attrib={
                                "attributeName":
                                setting["passthroughMapping"]["attributeName"],
                                "attributeType":
                                setting["passthroughMapping"]["attributeType"]
                            })
                    elif setting.get("discreteMapping"):
                        discrete_mapping = ET.SubElement(
                            visual_property,
                            "discreteMapping",
                            attrib={
                                "attributeName":
                                setting["discreteMapping"]
                                ["attributeName"].format(time=time),
                                "attributeType":
                                setting["discreteMapping"]["attributeType"]
                            })
                        for attribute_value, value in setting[
                                "discreteMapping"][
                                    "discreteMappingEntry"].items():
                            ET.SubElement(discrete_mapping,
                                          "discreteMappingEntry",
                                          attrib={
                                              "attributeValue":
                                              attribute_value,
                                              "value": value
                                          })

            for i, ptm in enumerate(
                    ppi_network.post_translational_modifications()[time]):
                for visual_property in visual_style.find("node").findall(
                        "visualProperty"):
                    if visual_property.get(
                            "name") == "NODE_CUSTOMGRAPHICS_{}".format(i + 1):
                        visual_property.set(
                            "default",
                            self.bar_chart(time,
                                           ptm,
                                           ppi_network.sites(time, ptm),
                                           cy_range=bar_chart_range))
                    elif visual_property.get(
                            "name"
                    ) == "NODE_CUSTOMGRAPHICS_POSITION_{}".format(i + 1):
                        visual_property.set("default",
                                            NODE_CUSTOMGRAPHICS_POSITION[i])

    def bar_chart(self, time, ptm, sites, cy_range=(-3.0, 3.0)):
        chart = json.dumps({
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
            ["{} {} {}".format(time, ptm, site + 1) for site in range(sites)]
        })
        return "org.cytoscape.BarChart: {}".format(chart)

    def export_as_xml(self, file_name):
        with open(file_name, "wb") as file:
            self.write(file, encoding="UTF-8", xml_declaration=True)
