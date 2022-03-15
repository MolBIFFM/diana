import math

LOGARITHM = {
    None:
        lambda changes, combination: math.log2(combination(changes)),
    2:
        lambda changes, combination: math.log2(
            combination([math.pow(change, 2.0) for change in changes])),
    10:
        lambda changes, combination: math.log2(
            combination([math.pow(change, 10.0) for change in changes])),
}
