EXTRACT = {
    ";":
    lambda entry: entry.split(";"),
    "|;":
    lambda entry:
    [e.split("|")[1] for e in entry.rstrip(";").split(";") if "|" in e],
    "|":
    lambda entry: [entry.split("|")[1]],
    "_,":
    lambda entry:
    [e.split("_")[0] for e in entry.rstrip(",").split(",") if "_" in e],
    "(":
    lambda entry: [entry.split("(")[0].rstrip()],
}
