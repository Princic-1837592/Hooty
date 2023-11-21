import os
import time
from argparse import ArgumentParser
from math import inf
from multiprocessing import cpu_count

import printers
from distances import K2P_distance, K2P_distance_ambiguity
from functions import read_species, read_fasta, compute_frequencies, remove_duplicates, compute_groups, run_processes, \
    compute_individual_groups


def main():
    start = time.time()

    parser = ArgumentParser()
    parser.add_argument("fasta_file")
    parser.add_argument("species_file")
    parser.add_argument(
        "-o",
        "--output",
        default="output.txt",
        metavar="output-file",
        dest="output_file"
    )
    parser.add_argument(
        "-l",
        "--log",
        default="log.txt",
        metavar="log-file",
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
    parser.add_argument(
        "--full-matrix",
        help="compute distances for all sequences",
        default=None,
        metavar="full-matrix",
        dest="full_matrix"
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

    dup_seqs = read_fasta(args.fasta_file, species, groups)
    print(f"Number of sequences with duplicates: {len(dup_seqs)}")

    frequencies = None
    distance_f = K2P_distance
    if args.ambiguous:
        frequencies = compute_frequencies(dup_seqs, n_groups, args.threshold)
        print("Computed frequencies")
        distance_f = K2P_distance_ambiguity

    seqs, seqs_offsets, min_dist_0 = remove_duplicates(dup_seqs, n_groups)
    print(f"Number of sequences without duplicates: {len(seqs)}")
    result = [[(inf, -inf)] * (i + 1) for i in range(n_groups)]
    processes = args.processes
    if processes > n_groups:
        processes = n_groups
    if processes == 1:
        p_result = compute_groups(
            range(n_groups),
            n_groups,
            None,
            min_dist_0,
            seqs,
            seqs_offsets,
            distance_f,
            frequencies
        )
        for g1 in range(n_groups):
            for g2_0, g2 in enumerate(range(g1, n_groups)):
                result[g2][g1] = p_result[g1][g2_0]
    else:
        run_processes(
            processes,
            n_groups,
            result,
            compute_groups,
            (min_dist_0, seqs, seqs_offsets, distance_f, frequencies)
        )

    string = printers.to_csv(result, species, groups)
    with open(args.output_file, "w") as f:
        f.write(string)
    print(f"Output written to {args.output_file}")

    if args.full_matrix is not None:
        n_groups = len(dup_seqs)
        processes = args.processes
        if processes > n_groups:
            processes = n_groups
        full_result = [[(inf, -inf)] * (i + 1) for i in range(n_groups)]
        run_processes(
            processes,
            n_groups,
            full_result,
            compute_individual_groups,
            (dup_seqs, distance_f, frequencies)
        )
        string = printers.to_csv(full_result, [seq.name for seq in dup_seqs], [i for i in range(n_groups)])
        with open(args.full_matrix, "w") as f:
            f.write(string)
        print(f"Full matrix written to {args.full_matrix}")

    end = time.time()
    print(f"{end - start:.2f}")


if __name__ == "__main__":
    main()
