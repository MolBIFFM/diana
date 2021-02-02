EXTRACT = {
    1: lambda entry: [entry],
    2: lambda entry: entry.split(";"),
    3: lambda entry: [entry.split("|")[1]],
    4: lambda entry: [
        subentry.split("_")[0]
        for subentry in entry.rstrip(",").split(",")
        if "_" in subentry
    ],
    5: lambda entry: [
        subentry.split("|")[1]
        for subentry in entry.rstrip(";").split(";")
        if "|" in subentry
    ],
}
