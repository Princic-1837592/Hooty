"""
    Hooty: fast DNA distance matrix
    Copyright (C) 2024  Andrea Princic & Giacomo Chiappa

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
"""

from dataclasses import dataclass
from multiprocessing import Process

from bases import REPLACEMENTS


@dataclass
class Sequence:
    index: int
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
