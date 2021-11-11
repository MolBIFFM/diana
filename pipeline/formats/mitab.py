def parse(entry):
    if entry == "-":
        return {}
    else:
        values = {}
        for namespace, identifier in (namespace_identifier.split(
                ":", 1) for namespace_identifier in entry.split("|")):

            if "(" in identifier and ")" in identifier:
                identifiers = (identifier[:identifier.find("(")].strip("\""),
                               identifier[identifier.find("(") +
                                          1:identifier.find(")")].strip("\""))
            else:
                identifiers = (identifier.strip("\""), None)

            if namespace in values:
                values[namespace].append(identifiers)
            else:
                values[namespace] = [identifiers]

        return values


def get_identifiers_from_namespace(entry, namespace):
    return [
        key for key, value in parse(entry).get(namespace, [("", None)]) if key
    ]


def get_terms_from_namespace(entry, namespace):
    return [
        value for key, value in parse(entry).get(namespace, [("", None)])
        if value
    ]


def namespace_has_identifier(entry, namespace, identifier):
    if values := parse(entry):
        return any(key == identifier for key, value in values[namespace])
    else:
        return False


def namespace_has_term(entry, namespace, term):
    if values := parse(entry):
        return any(value == term for key, value in values[namespace])
    else:
        return False


def namespace_has_any_identifier_from(entry, namespace, identifiers):
    if values := parse(entry):
        return any(key in identifiers for key, value in values[namespace])
    else:
        return False


def namespace_has_any_term_from(entry, namespace, terms):
    if values := parse(entry):
        return any(value in terms for key, value in values[namespace])
    else:
        return False
