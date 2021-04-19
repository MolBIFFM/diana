EXTRACT = {
    1: lambda entry: [entry],
    2: lambda entry: entry.split(";"),
    3: lambda entry: [
        subentry.split("|")[1]
        for subentry in entry.rstrip(";").split(";")
        if "|" in subentry
    ],
}
