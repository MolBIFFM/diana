MODIFY = {
    "X": lambda entry: entry,
    "X;Y": lambda entry: entry.split(";")[0],
    "Y|X|Y": lambda entry: entry.split("|")[1],
    "X_Y": lambda entry: entry.split("_")[0],
}
