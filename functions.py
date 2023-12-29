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

from math import inf
from multiprocessing import Pipe, Process

from bases import BYTES_MAP, REPLACEMENTS, T
from classes import AmbiguityInfo, Frequencies, Offset, ProcessData, Sequence


def read_species(species_file) -> tuple[list[str], list[int], list[int]]:
    """
    
    :param species_file: file containing the list groups of specie, one group per line, separated by commas
    :return: a tuple containing:
    - the list of all species
    - the list of groups, where the i-th element is the group of the i-th species
    - the list of offsets, where the j-th element is the index of the first species of the j-th group. The length of this list is the number of groups
    """
    species = []
    groups = []
    offsets = []
    with open(species_file, "r") as f:
        lines = f.read().strip().splitlines()
        group = 0
        for line in lines:
            offsets.append(len(species))
            for s in map(str.strip, line.split(",")):
                groups.append(group)
                species.append(s)
            group += 1
    return species, groups, offsets


def convert_to_single_lines(lines):
    buffer = []
    result = [lines[0]]
    lines.append(">")
    for l in lines[1:]:
        if l[0] == ">":
            result.append("".join(buffer))
            result.append(l)
            buffer.clear()
        else:
            buffer.append(l)
    return result


def to_bytes(sequence: str) -> bytes:
    return bytes(map(lambda c: BYTES_MAP[c], sequence))


def read_fasta(fasta_file, species, groups) -> list[Sequence]:
    with open(fasta_file, "r") as f:
        lines = f.read().strip().splitlines()
    # if len(lines) % 2 != 0:
    #     print("Error: fasta file is not in correct format")
    #     return None
    if lines[2][1] != ">":
        lines = convert_to_single_lines(lines)

    seqs = list()
    # fatal_error = False
    for l in range(0, len(lines), 2):
        name = lines[l][1:].strip()
        of_groups = set()
        for i, s in enumerate(species):
            if s in name:
                of_groups.add(groups[i])
        of_n_groups = len(of_groups)
        if of_n_groups == 1:
            of_group = of_groups.pop()
            seqs.append(Sequence(l, of_group, name, to_bytes(lines[l + 1].strip())))
        elif of_n_groups > 1:
            of_species = ", ".join(map(lambda g: species[g], of_groups))
            print(f"WARNING: species of sequence {name} is not unique. Found matches for: {of_species}")
            # fatal_error = True
        elif of_n_groups == 0 and len(name) > 0:
            print(f"WARNING: match not found for sequence '{name}'")
    # if fatal_error:
    #     exit(1)

    if any(len(seq.seq) != len(seqs[0].seq) for seq in seqs):
        print("ERROR: sequences in fasta file are not of equal length")
        exit(1)

    return seqs


def remove_duplicates(seqs: list[Sequence], n_groups: int) -> tuple[list[Sequence], list[Offset], list[bool]]:
    min_dist_0 = [False for _ in range(n_groups)]
    dedup = set()
    for seq in seqs:
        if seq not in dedup:
            dedup.add(seq)
        else:
            min_dist_0[seq.group] = True
    seqs = list(dedup)
    seqs.sort(key=lambda x: (x.group, x.index))
    offsets = [Offset(len(seqs), 0) for _ in range(n_groups)]
    for i, seq in enumerate(seqs):
        if offsets[seq.group].offset == len(seqs):
            offsets[seq.group].offset = i
        offsets[seq.group].count += 1
    return seqs, offsets, min_dist_0


def compute_frequencies(seqs: list[Sequence], n_groups, threshold) -> AmbiguityInfo:
    result = [[Frequencies() for _ in range(len(seqs[0].seq))] for _ in range(n_groups)]
    for seq in seqs:
        for i, base in enumerate(seq.seq):
            if base <= T:
                result[seq.group][i].count[base] += 1
                result[seq.group][i].normal_count += 1
            else:
                result[seq.group][i].ambiguous_count += 1
    for group in result:
        for freq in group:
            for i, replace in enumerate(REPLACEMENTS):
                rep_count = max(sum(freq.count[b] for b in replace), 1)
                for b in replace:
                    freq.frequencies[i][b] = freq.count[b] / rep_count
            freq.percentage = freq.ambiguous_count / max(freq.ambiguous_count + freq.normal_count, 1)
    return AmbiguityInfo(result, threshold)


def compute_groups(
        groups,
        n_groups: int,
        pipe,
        min_dist_0,
        seqs: list[Sequence],
        seqs_offsets: list[Offset],
        distance_f,
        frequencies: AmbiguityInfo
):
    result = [[(inf, -inf)] * (n_groups - group) for group in groups]
    for g1_0, g1 in enumerate(groups):
        g1_offset = seqs_offsets[g1]
        for g2_0, g2 in enumerate(range(g1, n_groups)):
            g2_offset = seqs_offsets[g2]
            min_score, max_score = inf, -inf
            for i in range(g1_offset.offset, g1_offset.offset + g1_offset.count):
                for j in range(g2_offset.offset, g2_offset.offset + g2_offset.count):
                    if i == j:
                        continue
                    distance = distance_f(seqs[i], seqs[j], frequencies)
                    if distance < min_score:
                        min_score = distance
                    if distance > max_score:
                        max_score = distance
            if g1 == g2:
                # fixme guardare qui per quel problema che mette male "0.0" e "/" sui gruppi
                if g1_offset.count == 1:
                    max_score = inf
                    min_score = inf
                elif min_dist_0[g1]:
                    min_score = 0.0
            result[g1_0][g2_0] = (min_score, max_score)
    if pipe is None:
        return result
    pipe.send(result)


def run_processes(n_processes: int, n_groups: int, result: list[list], target, args):
    numbers = list(range(n_groups))
    groups_per_chunk = n_groups // n_processes
    print(f"Working with {n_processes} processes, {groups_per_chunk} groups per process")
    chunks = [numbers[i:i + groups_per_chunk] for i in range(0, n_groups, groups_per_chunk)]
    if len(chunks) > n_processes:
        chunks[-2].extend(chunks[-1])
        chunks.pop()
    processes = []
    for g, chunk in enumerate(chunks):
        father, son = Pipe(False)
        p = Process(
            target=target,
            args=(chunk, n_groups, son, *args)
        )
        p.start()
        processes.append(ProcessData(p, father, son, chunk))
    for p, pdata in enumerate(processes):
        p_result = pdata.father.recv()
        for g1_0, g1 in enumerate(pdata.chunk):
            for g2_0, g2 in enumerate(range(g1, n_groups)):
                result[g2][g1] = p_result[g1_0][g2_0]
        pdata.father.close()
        pdata.son.close()
        pdata.process.join()


def compute_individual_groups(
        groups,
        n_groups: int,
        pipe,
        seqs: list[Sequence],
        distance_f,
        frequencies
):
    result = [[(inf, -inf)] * (n_groups - group) for group in groups]
    for g1_0, g1 in enumerate(groups):
        for g2_0, g2 in enumerate(range(g1, n_groups)):
            result[g1_0][g2_0] = distance_f(seqs[g1], seqs[g2], frequencies)
    if pipe is None:
        return result
    pipe.send(result)
