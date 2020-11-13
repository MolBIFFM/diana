def parse(entry):
    if entry != "-":
        values = {}
        for ns, identifier in (ns_identifier.split(":", 1)
                               for ns_identifier in entry.split("|")):
                                   
            identifier = identifier.split("(")
            identifier[0] = identifier[0].strip("\"")
            if len(identifier) == 2:
                identifier = identifier[1][:-1]

            if ns in values:
                values[ns].append(identifier)
            else:
                values[ns] = [identifier]
        return values
    else:
        return {}


def get_id_from_ns(entry, ns):
    return parse(entry).get(ns, [None])[0]


def ns_has_id(entry, ns, identifier):
    if values := parse(entry):
        return identifier in values[ns]
    else:
        return False


def ns_has_any_id(entry, ns, identifiers):
    if values := parse(entry):
        return any(identifier in values[ns] for identifier in identifiers)
    else:
        return False
