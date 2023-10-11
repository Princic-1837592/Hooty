import os
import sys
import time
from dataclasses import dataclass
from math import inf
from multiprocessing import cpu_count, Process, Pipe
from multiprocessing.connection import Connection
from typing import List, Tuple, Optional

import printers
from distances import K2Pdistance


@dataclass
class Sequence:
    group: int
    name: str
    seq: str

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
    father: Connection
    son: Connection
    chunk: List[int]


def read_species(species_file) -> Optional[Tuple[List[str], List[int], List[int]]]:
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


def read_fasta(fasta_file, species: List[str], groups: List[int]) -> Optional[
    Tuple[List[Sequence], List[Offset], List[bool]]
]:
    with open(fasta_file, "r") as f:
        lines = f.read().strip().splitlines()
    # if len(lines) % 2 != 0:
    #     print("Error: fasta file is not in correct format")
    #     return None
    if lines[2][1] != ">":
        lines = convert_to_single_lines(lines)

    seqs = set()
    min_dist_0 = [False for _ in range(groups[-1] + 1)]
    for l in range(0, len(lines), 2):
        name = lines[l][1:].strip()
        of_groups = set()
        for i, s in enumerate(species):
            if s in name:
                of_groups.add(groups[i])
        n_groups = len(of_groups)
        if n_groups == 1:
            old_len = len(seqs)
            of_group = of_groups.pop()
            seqs.add(Sequence(of_group, name, lines[l + 1].strip()))
            if len(seqs) == old_len:
                min_dist_0[of_group] = True
        elif n_groups > 1:
            pass
            # print("Warning: species of sequence is not unique")
        elif n_groups == 0:
            pass
            # print("Warning: species of sequence not found")
    seqs = list(seqs)

    if any(len(seq.seq) != len(seqs[0].seq) for seq in seqs):
        print("Error: sequences in fasta file are not of equal length")
        return None

    seqs.sort(key=lambda x: x.group)
    offsets = [Offset(len(seqs), 0) for _ in range(groups[-1] + 1)]
    for i, seq in enumerate(seqs):
        if offsets[seq.group].offset == len(seqs):
            offsets[seq.group].offset = i
        offsets[seq.group].count += 1
    return seqs, offsets, min_dist_0


def main():
    start = time.time()

    if len(sys.argv) < 3:
        print("Usage: python main.py <fasta_file> <species_file> [output_file] [processes]")
        return

    fasta_file, species_file = sys.argv[1], sys.argv[2]
    if not os.path.exists(fasta_file):
        print("Error: fasta file does not exist")
        return
    if not os.path.exists(species_file):
        print("Error: species file does not exist")
        return
    output_file = "output.txt"
    if len(sys.argv) > 3:
        output_file = sys.argv[3]
    else:
        print(f"Warning: no output file specified, using default {output_file}")
    p_count = cpu_count()
    if len(sys.argv) > 4:
        p_count = int(sys.argv[4])

    species, groups, group_offsets = read_species(species_file)
    n_groups = len(group_offsets)
    print(f"Number of species: {len(species)}")
    print(f"Number of groups: {n_groups}")
    if p_count > n_groups:
        p_count = n_groups

    seqs, seqs_offsets, min_dist_0 = read_fasta(fasta_file, species, groups)
    print(f"Number of sequences: {len(seqs)}")

    numbers = list(range(n_groups))
    groups_per_chunk = n_groups // p_count
    print(f"Working with {p_count} processes, {groups_per_chunk} groups per process")
    chunks = [numbers[i:i + groups_per_chunk] for i in range(0, n_groups, groups_per_chunk)]
    if len(chunks) > p_count:
        chunks[-2].extend(chunks[-1])
        chunks.pop()
    result = [[(inf, -inf)] * (i + 1) for i in range(n_groups)]
    processes = []
    for g, chunk in enumerate(chunks):
        father, son = Pipe(False)
        p = Process(target=compute_group, args=(chunk, n_groups, min_dist_0[g], seqs, seqs_offsets, son))
        p.start()
        processes.append(ProcessData(p, father, son, chunk))
    for p, pdata in enumerate(processes):
        print("receiving")
        p_result = pdata.father.recv()
        print("received", len(p_result), "results")
        for g1_0, g1 in enumerate(pdata.chunk):
            for g2_0, g2 in enumerate(range(g1, n_groups)):
                result[g2][g1] = p_result[g1_0][g2_0]
        pdata.father.close()
        pdata.son.close()
        pdata.process.join()
        print("joined")

    string = printers.to_csv(result, species, groups)
    with open(output_file, "w") as f:
        f.write(string)

    end = time.time()
    print(f"{end - start:.2f}")


def compute_group(groups, n_groups, min_dist_0, seqs, seqs_offsets, pipe):
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
                    distance = K2Pdistance(seqs[i].seq, seqs[j].seq)
                    if distance < min_score:
                        min_score = distance
                    if distance > max_score:
                        max_score = distance
            if g1 == g2 and min_dist_0:
                min_score = 0.0
            result[g1_0][g2_0] = (min_score, max_score)
    print(os.getpid(), "sending")
    pipe.send(result)
    print(os.getpid(), "sent")


if __name__ == "__main__":
    main()
