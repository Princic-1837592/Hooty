from math import nan, log, sqrt

import bases

TRANSITIONS = [
    # A  C  G  T
    [0, 0, 1, 0],
    [0, 0, 0, 1],
    [1, 0, 0, 0],
    [0, 1, 0, 0],
]


# https://github.com/kgori/python_tools_on_github/blob/master/pairwise_distances.py
def K2P_distance(seq1, seq2, _):
    """
    Kimura 2-Parameter distance = -0.5 log( (1 - 2p -q) * sqrt( 1 - 2q ) )
    where:
    p = transition frequency
    q = transversion frequency
    """
    pairs = []

    # collect ungapped pairs
    for x in zip(seq1, seq2):
        if bases.GAP not in x and bases.N not in x:
            pairs.append(x)

    ts_count = 0
    tv_count = 0
    length = len(pairs)

    for (x, y) in filter(lambda p: p[0] != p[1] and p[0] <= bases.T and p[1] <= bases.T, pairs):
        ts_count += TRANSITIONS[x][y]
        tv_count += 1 - TRANSITIONS[x][y]

    p = float(ts_count) / length
    q = float(tv_count) / length
    try:
        d = -0.5 * log((1 - 2 * p - q) * sqrt(1 - 2 * q))
    except ValueError:
        print("Tried to take log of a negative number")
        return nan
    return d


def K2P_distance_ambiguity(seq1, seq2, frequencies):
    """
    Kimura 2-Parameter distance = -0.5 log( (1 - 2p -q) * sqrt( 1 - 2q ) )
    where:
    p = transition frequency
    q = transversion frequency
    """
    pairs = []

    # collect ungapped pairs
    for x in zip(seq1, seq2):
        if bases.GAP not in x and bases.N not in x:
            pairs.append(x)

    ts_count = 0
    tv_count = 0
    length = len(pairs)

    for (x, y) in filter(lambda p: p[0] != p[1] and p[0] <= bases.T and p[1] <= bases.T, pairs):
        ts_count += TRANSITIONS[x][y]
        tv_count += 1 - TRANSITIONS[x][y]

    p = float(ts_count) / length
    q = float(tv_count) / length
    try:
        d = -0.5 * log((1 - 2 * p - q) * sqrt(1 - 2 * q))
    except ValueError:
        print("Tried to take log of a negative number")
        return nan
    return d
