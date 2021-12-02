def parse(entry):
    if entry == "-":
        return {}
    else:
        values = {}
        for namespace, identifier_term in (namespace_identifier.split(
                ":", 1) for namespace_identifier in entry.split("|")):

            if "(" in identifier_term and ")" in identifier_term:
                identifier = identifier_term[:identifier_term.find("(")].strip(
                    "\"")

                term = identifier_term[identifier_term.find("(") +
                                       1:identifier_term.find(")")].strip("\"")
            else:
                identifier = identifier_term.strip("\"")
                term = None

            if namespace in values:
                values[namespace].append((identifier, term))
            else:
                values[namespace] = [(identifier, term)]

        return values


def get_identifiers_from_namespace(entry, namespace):
    return [
        identifier
        for identifier, _ in parse(entry).get(namespace, [("", None)])
        if identifier
    ]


def get_terms_from_namespace(entry, namespace):
    return [
        term for _, term in parse(entry).get(namespace, [("", None)]) if term
    ]


def namespace_has_identifier(entry, namespace, search_identifier):
    if values := parse(entry):
        return any(search_identifier == identifier
                   for identifier, _ in values[namespace])
    else:
        return False


def namespace_has_term(entry, namespace, search_term):
    if values := parse(entry):
        return any(search_term == term for _, term in values[namespace])
    else:
        return False


def namespace_has_any_identifier_from(entry, namespace, identifiers):
    if values := parse(entry):
        return any(identifier in identifiers
                   for identifier, _ in values[namespace])
    else:
        return False


def namespace_has_any_term_from(entry, namespace, terms):
    if values := parse(entry):
        return any(term in terms for _, term in values[namespace])
    else:
        return False
