def parse(entry):
    if entry != "-":
        values = {}
        for ns, identifier in (
            ns_identifier.split(":", 1) for ns_identifier in entry.split("|")
        ):

            identifiers = identifier.split("(")
            identifiers = {
                identifiers[0].strip('"'): identifiers[1].rstrip(")")
                if len(identifiers) == 2
                else None
            }

            if ns in values:
                values[ns].append(identifiers)
            else:
                values[ns] = [identifiers]

        return values
    else:
        return {}


def get_id_from_ns(entry, ns):
    return list(parse(entry).get(ns, [{None: ""}])[0].keys())[0]


def ns_has_id(entry, ns, identifier):
    if values := parse(entry):
        return any(identifier in value.keys() for value in values[ns])
    else:
        return False


def ns_has_term(entry, ns, identifier):
    if values := parse(entry):
        return any(identifier in value.values() for value in values[ns])
    else:
        return False


def ns_has_any_id(entry, ns, identifiers):
    if values := parse(entry):
        return any(
            identifier in value.keys()
            for value in values[ns]
            for identifier in identifiers
        )
    else:
        return False


def ns_has_any_term(entry, ns, identifiers):
    if values := parse(entry):
        return any(
            identifier in value.values()
            for value in values[ns]
            for identifier in identifiers
        )
    else:
        return False
