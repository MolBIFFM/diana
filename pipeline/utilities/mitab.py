def parse(entry):
    if entry == "-":
        return {}
    else:
        values = {}
        for namespace, identifier in (namespace_identifier.split(
                ":", 1) for namespace_identifier in entry.split("|")):

            identifiers = identifier.split("(")
            identifiers = {
                identifiers[0].strip('"'):
                identifiers[1].rstrip(")") if len(identifiers) == 2 else None
            }

            if namespace in values:
                values[namespace].append(identifiers)
            else:
                values[namespace] = [identifiers]

        return values


def get_identifier_from_namespace(entry, namespace):
    return list(parse(entry).get(namespace, [{None: ""}])[0].keys())[0]


def namespace_has_identifier(entry, namespace, identifier):
    if values := parse(entry):
        return any(
            str(identifier) in value.keys() for value in values[namespace])
    else:
        return False


def namespace_has_term(entry, namespace, term):
    if values := parse(entry):
        return any(term in value.values() for value in values[namespace])
    else:
        return False


def namespace_has_any_identifier_from(entry, namespace, identifiers):
    if values := parse(entry):
        return any(
            str(identifier) in value.keys() for value in values[namespace]
            for identifier in identifiers)
    else:
        return False


def namespace_has_any_term_from(entry, namespace, terms):
    if values := parse(entry):
        return any(term in value.values() for value in values[namespace]
                   for term in terms)
    else:
        return False
