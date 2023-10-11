import os
import sys
import time
from dataclasses import dataclass
from math import inf
from typing import List, Tuple, Optional

from distances import K2Pdistance
import printers


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
        print("Usage: python main.py <fasta_file> <species_file> [output_file]")
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

    species, groups, group_offsets = read_species(species_file)
    n_groups = len(group_offsets)
    print(f"Number of species: {len(species)}")
    print(f"Number of groups: {n_groups}")

    seqs, seqs_offsets, min_dist_0 = read_fasta(fasta_file, species, groups)
    print(f"Number of sequences: {len(seqs)}")

    result = [[(inf, -inf)] * (i + 1) for i, _ in enumerate(range(len(group_offsets)))]
    for g1 in range(len(group_offsets)):
        group = g1 + 1
        perc = group / n_groups * 100
        eta = (time.time() - start) / group * (n_groups - g1 - 1)
        print(f"Group {group} ({perc:.2f}% ETA: {eta:.2f}s)")
        for g2 in range(g1, len(group_offsets)):
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
            if g1 == g2 and min_dist_0[g1]:
                min_score = 0.0
            result[g2][g1] = (min_score, max_score)

    string = printers.to_csv(result, species, groups)
    with open(output_file, "w") as f:
        f.write(string)

    end = time.time()
    print(f"{end - start:.2f}")


if __name__ == "__main__":
    main()
