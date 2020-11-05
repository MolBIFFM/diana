def parse(entry):
    if entry != "-":
        values = {}
        for ns, identifiers in (ns_identifiers.split(":", 1)
                                for ns_identifiers in entry.split("|")):
            identifiers = identifiers.split("(")
            identifiers[0] = identifiers[0].strip("\"")
            if len(identifiers) == 2:
                identifiers[1] = identifiers[1][:-1]

            if ns in values:
                values[ns].extend([
                    identifier for identifier in identifiers
                    if identifier not in values[ns]
                ])
            else:
                values[ns] = identifiers
        return values
    else:
        return {}


def get_id_from_ns(entry, ns):
    return parse(entry).get(ns, [None])[0] if entry != "-" else None


def ns_has_id(entry, ns, identifier):
    if (values := parse(entry)):
        return identifier in values[ns]
    else:
        return False


def ns_has_any_id(entry, ns, identifiers):
    if (values := parse(entry)):
        return any(identifier in values[ns] for identifier in identifiers)
    else:
        return False
