EXTRACT = {
    "X": lambda entry: [entry],
    "X;X": lambda entry: entry.split(";"),
    "Y|X|Z": lambda entry: [entry.split("|")[1]],
    "X_Y,X_Y": lambda entry: [
        subentry.split("_")[0] for subentry in entry.split(",") if "_" in subentry
    ],
    "Y|X|Z;Y|X|Z": lambda entry: [
        subentry.split("|")[1]
        for subentry in entry.rstrip(";").split(";")
        if "|" in subentry
    ],
}
