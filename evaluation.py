import os
import time
from argparse import ArgumentParser
from math import inf
from multiprocessing import Pipe, Process, cpu_count

from classes import ProcessData
from distances import K2P_distance, K2P_distance_ambiguity
from main import compute_frequencies, compute_groups, read_fasta, read_species, remove_duplicates
from printers import format_value


def to_columns(matrix, min_max, dataset) -> list[str]:
    min_max = 0 if min_max == "min" else 1
    lines = [dataset]
    for row in matrix:
        for vals in row:
            lines.append(format_value(vals[min_max]))
    return lines


def main():
    start = time.time()

    parser = ArgumentParser()
    parser.add_argument("directory")
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
        dest="processes",
        type=int
    )
    parser.add_argument(
        "-a",
        "--ambiguous",
        default=False,
        action="store_true",
        help="count ambiguous bases",
        dest="ambiguous"
    )
    parser.add_argument(
        "-t",
        "--threshold",
        default=0.1,
        help="max percentage of ambiguous bases in a group",
        dest="threshold",
        type=float
    )
    args = parser.parse_args()

    fastas = []
    for file in os.listdir(args.directory):
        if file.endswith(".fasta"):
            fastas.append(os.path.join(args.directory, file))
    fastas.sort()
    species, groups, group_offsets = read_species(args.species_file)
    n_groups = len(group_offsets)
    print(f"Number of species: {len(species)}")
    print(f"Number of groups: {n_groups}")
    if args.processes > n_groups:
        args.processes = n_groups

    columns_min = []
    columns_max = []
    for file in fastas:
        print(f"Reading file: {file}")
        seqs = read_fasta(file, species, groups)
        print(f"Number of sequences with duplicates: {len(seqs)}")

        frequencies = None
        distance_f = K2P_distance
        if args.ambiguous:
            frequencies = compute_frequencies(seqs, n_groups, args.threshold)
            print("Computed frequencies")
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

        columns_min.append(to_columns(result, "min", os.path.splitext(os.path.basename(file))[0]))
        columns_max.append(to_columns(result, "max", os.path.splitext(os.path.basename(file))[0]))
        print("\n")
    name, ext = os.path.splitext(args.output_file)
    with open(f"{name}_min{ext}", "w") as f:
        for i in range(len(columns_min[0])):
            f.write("\t".join(c[i] for c in columns_min))
            f.write("\n")
    with open(f"{name}_max{ext}", "w") as f:
        for i in range(len(columns_max[0])):
            f.write("\t".join(c[i] for c in columns_max))
            f.write("\n")

    end = time.time()
    print(f"{end - start:.2f}")


if __name__ == "__main__":
    main()
