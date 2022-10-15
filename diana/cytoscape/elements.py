"""Additions to Cytoscape style element trees."""
import xml.etree.ElementTree as ET
from typing import Literal


def add_dependency(parent: ET.Element, name: str, value: bool) -> ET.Element:
    """
    Attaches a dependency to a parent element of a Cytoscape style.

    Args:
        parent: The parent element to attach the dependency to.
        name: The name of the dependency.
        value: The value of the dependency.

    Returns:
        The dependency element.
    """
    return ET.SubElement(parent,
                         "dependency",
                         attrib={
                             "name": name,
                             "value": {
                                 False: "false",
                                 True: "true"
                             }[value]
                         })


def add_visual_property(parent: ET.Element, name: str,
                        default: bool | int | str | float) -> ET.Element:
    """
    Attaches a visual property to a parent element of a Cytoscape style.

    Args:
        parent: The parent element to attach the visual property to.
        name: The name of the visual property.
        default: The default value of the visual property.

    Returns:
        The visual property element.
    """
    return ET.SubElement(
        parent,
        "visualProperty",
        attrib={
            "name": name,
            "default": {
                False: "false",
                True: "true"
            }[default] if isinstance(default, bool) else str(default)
        })


def add_continuous_mapping(
    parent: ET.Element, attribute_name: str, attribute_type: Literal["float",
                                                                     "integer",
                                                                     "string"],
    entries: dict[int | float, tuple[int | str | float, int | str | float,
                                     int | str | float]]
) -> ET.Element:
    """
    Attaches a continuous mapping to a visual property element of a Cytoscape
    style.

    Args:
        parent: The visual property element to attach the continuous mapping to.
        attribute_name: The name of the network attribute the visual property
            depends on.
        default: The type of the network attribute the visual property depends
            on.
        entries: A mapping from continuous values of the network attribute to
            the value of the visual property.

    Returns:
        The visual property element the continuous mapping is attached to.
    """
    continuous_mapping = ET.SubElement(parent,
                                       "continuousMapping",
                                       attrib={
                                           "attributeName": attribute_name,
                                           "attributeType": attribute_type
                                       })
    for key, (lesser_value, equal_value, greater_value) in entries.items():
        ET.SubElement(continuous_mapping,
                      "continuousMappingPoint",
                      attrib={
                          "attrValue": str(key),
                          "lesserValue": str(lesser_value),
                          "equalValue": str(equal_value),
                          "greaterValue": str(greater_value)
                      })
    return continuous_mapping


def add_discrete_mapping(parent: ET.Element, attribute_name: str,
                         attribute_type: Literal["float", "integer", "string"],
                         entries: dict[str, str]) -> ET.Element:
    """
    Attaches a discrete mapping to a visual property element of a Cytoscape
    style.

    Args:
        parent: The visual property element to attach the discrete mapping to.
        attribute_name: The name of the network attribute the visual property
            depends on.
        default: The type of the network attribute the visual property depends
            on.
        entries: A mapping from discrete values of the network attribute to
            the value of the visual property.

    Returns:
        The visual property element the discrete mapping is attached to.
    """
    discrete_mapping = ET.SubElement(parent,
                                     "discreteMapping",
                                     attrib={
                                         "attributeName": attribute_name,
                                         "attributeType": attribute_type
                                     })

    for key, value in entries.items():
        ET.SubElement(discrete_mapping,
                      "discreteMappingEntry",
                      attrib={
                          "attributeValue": str(key),
                          "value": str(value)
                      })

    return discrete_mapping


def add_passthrough_mapping(
        parent: ET.Element, attribute_name: str,
        attribute_type: Literal["float", "integer", "string"]) -> ET.Element:
    """
    Attaches a passthrough mapping to a visual property element of a Cytoscape
    style.

    Args:
        parent: The visual property element to attach the passthrough mapping
            to.
        attribute_name: The name of the network attribute the visual property
            depends on.
        default: The type of the network attribute the visual property depends
            on.

    Returns:
        The visual property element the passthrough mapping is attached to.
    """
    return ET.SubElement(parent,
                         "passthroughMapping",
                         attrib={
                             "attributeName": attribute_name,
                             "attributeType": attribute_type
                         })
