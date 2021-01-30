EXTRACT = {
    "X": lambda entry: [entry],
    "X;X": lambda entry: entry.split(";"),
    "Y|X|Z": lambda entry: entry.split("|")[1:2],
    "X_Y,X_Y": lambda entry: entry.split(",")[0].split("_")[:1],
    "Y|X|Z;Y|X|Z": lambda entry: [
        subentry.split("|")[1]
        for subentry in entry.rstrip(";").split(";")
        if len(subentry.split("|")) > 1
    ],
}
