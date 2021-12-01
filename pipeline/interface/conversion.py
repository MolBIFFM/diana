import math

LOGARITHM = {
    None: math.log2,
    2: lambda x: x,
    10: lambda x: math.log10(x) / math.log10(2.0),
}
