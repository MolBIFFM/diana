MODIFY = {
    "X": lambda x: x,
    "X.Y;": lambda x: x.split(";")[0].split(".")[0],
    "Y|X|Y": lambda x: x.split("|")[1],
    "X_Y": lambda x: x.split("_")[0],
}
