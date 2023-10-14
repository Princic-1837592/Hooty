from math import log, nan, sqrt

from bases import A, C, G, GAP, N, T, W
from classes import Frequencies, Sequence

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

    p = float(ts_count) / length
    q = float(tv_count) / length
    try:
        d = -0.5 * log((1 - 2 * p - q) * sqrt(1 - 2 * q))
    except ValueError:
        print("Tried to take log of a negative number")
        return nan
    return d


def K2P_distance_ambiguity(seq1: Sequence, seq2: Sequence, frequencies: list[list[Frequencies]]):
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
    f1 = frequencies[seq1.group]
    f2 = frequencies[seq2.group]

    for i, (x, y) in filter(lambda p: p[1][0] != p[1][1], pairs):
        if x <= T and y <= T:
            ts_count += TRANSITIONS[x][y]
            tv_count += 1 - TRANSITIONS[x][y]
        else:
            x -= W
            y -= W
            ts_count_f += (f1[i].frequencies[x][A] * f2[i].frequencies[y][G] +
                           f1[i].frequencies[x][G] * f2[i].frequencies[y][A] +
                           f1[i].frequencies[x][C] * f2[i].frequencies[y][T] +
                           f1[i].frequencies[x][T] * f2[i].frequencies[y][C])
            tv_count_f += (f1[i].frequencies[x][A] * f2[i].frequencies[y][T] +
                           f1[i].frequencies[x][T] * f2[i].frequencies[y][A] +
                           f1[i].frequencies[x][A] * f2[i].frequencies[y][C] +
                           f1[i].frequencies[x][C] * f2[i].frequencies[y][A] +
                           f1[i].frequencies[x][G] * f2[i].frequencies[y][T] +
                           f1[i].frequencies[x][T] * f2[i].frequencies[y][G] +
                           f1[i].frequencies[x][G] * f2[i].frequencies[y][C] +
                           f1[i].frequencies[x][C] * f2[i].frequencies[y][G])
    p = (float(ts_count) + ts_count_f) / length
    q = (float(tv_count) + tv_count_f) / length
    try:
        d = -0.5 * log((1 - 2 * p - q) * sqrt(1 - 2 * q))
    except ValueError:
        print("Tried to take log of a negative number")
        return nan
    return d
