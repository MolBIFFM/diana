import xml.etree.ElementTree as ET

import pipeline.configuration.cytoscape_style


class CytoscapeStyle(ET.ElementTree):
    def __init__(self, times):
        super(CytoscapeStyle, self).__init__(
            ET.Element("vizmap",
                       attrib={
                           "id": "VizMap",
                           "documentVersion": "3.0"
                       }))
        for time in times:
            visual_style = ET.SubElement(self.getroot(),
                                         "visualStyle",
                                         attrib={"name": str(time)})
            tags = {}
            for tag in pipeline.configuration.cytoscape_style.settings:
                tags[tag] = ET.SubElement(visual_style, tag)
                for name, setting in pipeline.configuration.cytoscape_style.settings[
                        tag]["dependency"].items():
                    ET.SubElement(tags[tag],
                                  "dependency",
                                  attrib={
                                      "name": name,
                                      "value": setting["value"]
                                  })
                for name, setting in pipeline.configuration.cytoscape_style.settings[
                        tag]["visualProperty"].items():
                    subtag = ET.SubElement(tags[tag],
                                           "visualProperty",
                                           attrib={
                                               "name": name,
                                               "default": setting["default"]
                                           })
                    if setting.get("passthroughMapping"):
                        ET.SubElement(
                            subtag,
                            "passthroughMapping",
                            attrib={
                                "attributeName":
                                setting["passthroughMapping"]["attributeName"],
                                "attributeType":
                                setting["passthroughMapping"]["attributeType"]
                            })

    def export_as_xml(self, file_name):
        with open(file_name, "wb") as file:
            self.write(file, encoding="UTF-8", xml_declaration=True)
