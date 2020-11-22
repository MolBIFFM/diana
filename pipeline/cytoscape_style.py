import xml.etree.ElementTree as ET

import pipeline.configuration.cytoscape_style


class CytoscapeStyle(ET.ElementTree):
    def __init__(self, ppi_network):
        super(CytoscapeStyle, self).__init__(
            ET.Element("vizmap",
                       attrib={
                           "id": "VizMap",
                           "documentVersion": "3.0"
                       }))
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

    def export_as_xml(self, file_name):
        with open(file_name, "wb") as file:
            self.write(file, encoding="UTF-8", xml_declaration=True)
