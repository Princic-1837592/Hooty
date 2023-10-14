import os
import time
from argparse import ArgumentParser
from dataclasses import dataclass
from math import inf
from multiprocessing import Pipe, Process, cpu_count

import bases
import printers
from bases import BYTES_MAP
from distances import K2P_distance, K2P_distance_ambiguity


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
    frequencies: list[float]
    normal_count: int
    ambiguous_count: int

    def __init__(self):
        self.count = [0 for _ in range(4)]
        self.frequencies = []
        self.normal_count = 0
        self.ambiguous_count = 0


def read_species(species_file) -> tuple[list[str], list[int], list[int]] | None:
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
    for l in range(0, len(lines), 2):
        name = lines[l][1:].strip()
        of_groups = set()
        for i, s in enumerate(species):
            if s in name:
                of_groups.add(groups[i])
        of_n_groups = len(of_groups)
        if of_n_groups == 1:
            of_group = of_groups.pop()
            seqs.append(Sequence(of_group, name, to_bytes(lines[l + 1].strip())))
        elif of_n_groups > 1:
            of_species = ", ".join(map(lambda g: species[g], of_groups))
            print(f"Error: species of sequence is not unique. Found matches for: {of_species}")
            # todo
            # exit(1)
        elif of_n_groups == 0:
            # print(f"Warning: match not found for sequence '{name}'")
            pass

    if any(len(seq.seq) != len(seqs[0].seq) for seq in seqs):
        print("Error: sequences in fasta file are not of equal length")
        exit(1)

    return seqs


def remove_duplicates(seqs, n_groups) -> tuple[list[Sequence], list[Offset], list[bool]]:
    min_dist_0 = [False for _ in range(n_groups)]
    dedup = set()
    for seq in seqs:
        if seq not in dedup:
            dedup.add(seq)
        else:
            min_dist_0[seq.group] = True
    seqs = list(dedup)
    seqs.sort(key=lambda x: x.group)
    offsets = [Offset(len(seqs), 0) for _ in range(n_groups)]
    for i, seq in enumerate(seqs):
        if offsets[seq.group].offset == len(seqs):
            offsets[seq.group].offset = i
        offsets[seq.group].count += 1
    return seqs, offsets, min_dist_0


def compute_frequencies(seqs: list[Sequence], n_groups) -> list[list[Frequencies]]:
    result = [[Frequencies() for _ in range(len(seqs[0].seq))] for _ in range(n_groups)]
    for seq in seqs:
        for i, base in enumerate(seq.seq):
            if base <= bases.T:
                result[seq.group][i].count[base] += 1
                result[seq.group][i].normal_count += 1
            elif base < bases.N:
                result[seq.group][i].ambiguous_count += 1
    for group in result:
        for freq in group:
            d = max(freq.normal_count, 1)
            freq.frequencies = [c / d for c in freq.count]
    return result


def main():
    start = time.time()

    parser = ArgumentParser()
    parser.add_argument("fasta_file")
    parser.add_argument("species_file")
    parser.add_argument(
        "-o",
        "--output",
        default="output.txt",
        metavar="output file",
        dest="output_file"
    )
    parser.add_argument(
        "-l",
        "--log",
        default="log.txt",
        metavar="log file",
        dest="log_file"
    )
    parser.add_argument(
        "-p",
        "--processes",
        default=cpu_count(),
        metavar="processes",
        dest="processes"
    )
    parser.add_argument(
        "-a",
        "--ambiguous",
        default=False,
        action="store_true",
        help="count ambiguous bases",
        dest="ambiguous"
    )
    args = parser.parse_args()

    if not os.path.exists(args.fasta_file):
        print("Error: fasta file does not exist")
        return
    if not os.path.exists(args.species_file):
        print("Error: species file does not exist")
        return

    species, groups, group_offsets = read_species(args.species_file)
    n_groups = len(group_offsets)
    print(f"Number of species: {len(species)}")
    print(f"Number of groups: {n_groups}")
    if args.processes > n_groups:
        args.processes = n_groups

    seqs = read_fasta(args.fasta_file, species, groups)
    print(f"Number of sequences with duplicates: {len(seqs)}")

    frequencies = None
    distance_f = K2P_distance
    if args.ambiguous:
        frequencies = compute_frequencies(seqs, n_groups)
        distance_f = K2P_distance_ambiguity

    seqs, seqs_offsets, min_dist_0 = remove_duplicates(seqs, n_groups)
    print(f"Number of sequences without duplicates: {len(seqs)}")
    result = [[(inf, -inf)] * (i + 1) for i in range(n_groups)]
    if args.processes == 1:
        p_result = compute_groups(
            range(n_groups),
            n_groups,
            min_dist_0,
            seqs,
            seqs_offsets,
            None,
            distance_f,
            frequencies
        )
        for g1 in range(n_groups):
            for g2_0, g2 in enumerate(range(g1, n_groups)):
                result[g2][g1] = p_result[g1][g2_0]
    else:
        numbers = list(range(n_groups))
        groups_per_chunk = n_groups // args.processes
        print(f"Working with {args.processes} processes, {groups_per_chunk} groups per process")
        chunks = [numbers[i:i + groups_per_chunk] for i in range(0, n_groups, groups_per_chunk)]
        if len(chunks) > args.processes:
            chunks[-2].extend(chunks[-1])
            chunks.pop()
        processes = []
        for g, chunk in enumerate(chunks):
            father, son = Pipe(False)
            p = Process(
                target=compute_groups,
                args=(chunk, n_groups, min_dist_0, seqs, seqs_offsets, son, distance_f, frequencies)
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

    string = printers.to_csv(result, species, groups)
    with open(args.output_file, "w") as f:
        f.write(string)

    end = time.time()
    print(f"{end - start:.2f}")


def compute_groups(groups, n_groups, min_dist_0, seqs, seqs_offsets, pipe, distance_f, frequencies):
    result = [[(inf, -inf)] * (n_groups - group) for group in groups]
    for g1_0, g1 in enumerate(groups):
        for g2_0, g2 in enumerate(range(g1, n_groups)):
            g1_offset = seqs_offsets[g1]
            g2_offset = seqs_offsets[g2]
            min_score, max_score = inf, -inf
            for i in range(g1_offset.offset, g1_offset.offset + g1_offset.count):
                for j in range(g2_offset.offset, g2_offset.offset + g2_offset.count):
                    if i == j:
                        continue
                    distance = distance_f(seqs[i].seq, seqs[j].seq, frequencies)
                    if distance < min_score:
                        min_score = distance
                    if distance > max_score:
                        max_score = distance
            if g1 == g2:
                if min_dist_0[g1]:
                    min_score = 0.0
                if g1_offset.count == 1:
                    max_score = 0.0
            result[g1_0][g2_0] = (min_score, max_score)
    if pipe is None:
        return result
    pipe.send(result)


if __name__ == "__main__":
    main()
