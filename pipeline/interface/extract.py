EXTRACT = {
    1: lambda entry: [entry],
    2: lambda entry: entry.split(";"),
    3: lambda entry: [
        e.split("|")[1] for e in entry.rstrip(";").split(";") if "|" in e
    ],
    4: lambda entry: [entry.split("|")[1]],
    5: lambda entry: [
        e.split("_")[0] for e in entry.rstrip(",").split(",") if "_" in e
    ],
    6: lambda entry: [entry.split("(")[0].rstrip()],
}
