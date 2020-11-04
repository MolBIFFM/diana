def parse(entry):
    if entry != "-":
        values = {}
        for ns, ids in (ns_ids.split(":", 1) for ns_ids in entry.split("|")):
            ids = ids.split("(")
            ids[0] = ids[0].strip("\"")
            if len(ids) == 2:
                ids[1] = ids[1][:-1]
            
            if ns in values:
                values[ns].extend(ids)
            else:
                values[ns] = ids
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