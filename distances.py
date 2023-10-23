from math import log, nan, sqrt

from bases import A, C, G, GAP, N, T, W
from classes import AmbiguityInfo, Sequence

TRANSITIONS = [
    # A  C  G  T
    [0, 0, 1, 0],
    [0, 0, 0, 1],
    [1, 0, 0, 0],
    [0, 1, 0, 0],
]


# https://github.com/kgori/python_tools_on_github/blob/master/pairwise_distances.py
def K2P_distance(seq1: Sequence, seq2: Sequence, _):
    """
    Kimura 2-Parameter distance = -0.5 log( (1 - 2p -q) * sqrt( 1 - 2q ) )
    where:
    p = transition frequency
    q = transversion frequency
    """
    pairs = []

    # collect ungapped pairs
    for x in zip(seq1.seq, seq2.seq):
        if GAP not in x and N not in x:
            pairs.append(x)

    ts_count = 0
    tv_count = 0
    length = len(pairs)

    for (x, y) in filter(lambda p: p[0] != p[1] and p[0] <= T and p[1] <= T, pairs):
        ts_count += TRANSITIONS[x][y]
        tv_count += 1 - TRANSITIONS[x][y]

    p = float(ts_count) / max(length, 1)
    q = float(tv_count) / max(length, 1)
    try:
        d = -0.5 * log((1 - 2 * p - q) * sqrt(1 - 2 * q))
    except ValueError:
        print("Tried to take log of a negative number")
        return nan
    return d


def K2P_distance_ambiguity(seq1: Sequence, seq2: Sequence, frequencies: AmbiguityInfo):
    """
    Kimura 2-Parameter distance = -0.5 log( (1 - 2p -q) * sqrt( 1 - 2q ) )
    Kimura 2-Parameter distance = -0.5 log( (1 - 2p -q) * sqrt( 1 - 2q ) )
    where:
    p = transition frequency
    q = transversion frequency
    """
    pairs = []

    # collect ungapped pairs
    for i, x in enumerate(zip(seq1.seq, seq2.seq)):
        if GAP not in x and N not in x:
            pairs.append((i, x))

    ts_count = 0
    ts_count_f = 0.0
    tv_count = 0
    tv_count_f = 0.0
    length = len(pairs)
    g1 = frequencies[seq1.group]
    g2 = frequencies[seq2.group]

    for i, (x, y) in filter(lambda p: p[1][0] != p[1][1], pairs):
        if x <= T and y <= T:
            ts_count += TRANSITIONS[x][y]
            tv_count += 1 - TRANSITIONS[x][y]
        elif g1[i].percentage < frequencies.threshold and g2[i].percentage < frequencies.threshold:
            x = g1[i][x - W]
            y = g2[i][y - W]
            ts_count_f += x[A] * y[G] + x[G] * y[A] + x[C] * y[T] + x[T] * y[C]
            tv_count_f += (x[A] * y[T] + x[T] * y[A] +
                           x[A] * y[C] + x[C] * y[A] +
                           x[G] * y[T] + x[T] * y[G] +
                           x[G] * y[C] + x[C] * y[G])
        else:
            length -= 1
    p = (float(ts_count) + ts_count_f) / max(length, 1)
    q = (float(tv_count) + tv_count_f) / max(length, 1)
    try:
        d = -0.5 * log((1 - 2 * p - q) * sqrt(1 - 2 * q))
    except ValueError:
        print("Tried to take log of a negative number")
        return nan
    return d
