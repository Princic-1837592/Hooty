from dataclasses import dataclass
from multiprocessing import Process

from bases import REPLACEMENTS


@dataclass
class Sequence:
    group: int
    name: str
    seq: bytes

    def __hash__(self):
        return hash((self.group, self.seq))

    def __eq__(self, other):
        return self.group == other.group and self.seq == other.seq


@dataclass
class Offset:
    offset: int
    count: int


@dataclass
class ProcessData:
    process: Process
    father: ...
    son: ...
    chunk: list[int]


@dataclass
class Frequencies:
    count: list[int]
    frequencies: list[list[float]]
    normal_count: int
    ambiguous_count: int
    percentage: float

    def __init__(self):
        self.count = [0 for _ in range(4)]
        self.frequencies = [[0.0, 0.0, 0.0, 0.0] for _ in range(len(REPLACEMENTS))]
        self.normal_count = 0
        self.ambiguous_count = 0
        self.percentage = 0.0

    def __getitem__(self, item):
        return self.frequencies[item]


@dataclass
class AmbiguityInfo:
    groups: list[list[Frequencies]]
    threshold: float

    def __getitem__(self, item):
        return self.groups[item]
