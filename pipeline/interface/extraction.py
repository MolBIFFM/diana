import itertools

EXTRACT = {
    ";":
    lambda entry: entry.split(";"),
    "|;":
    lambda entry: [
        site.split("|")[1] for site in entry.rstrip(";").split(";")
        if "|" in site
    ],
    "|":
    lambda entry: [entry.split("|")[1]],
    "_,":
    lambda entry: [
        site.split("_")[0] for site in entry.rstrip(",").split(",")
        if "_" in site
    ],
    "[();()]":
    lambda entry: [
        position for position in itertools.chain.from_iterable([
            int(site[1:site.find("(")]) if "(" in site and ")" else None
            for site in subentry[subentry.find("[") + 1:subentry.find("]")].
            split(";")
        ] for subentry in entry.split(";")) if position is not None
    ]
}
