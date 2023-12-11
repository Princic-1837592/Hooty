import os
import time
from math import inf
from multiprocessing import cpu_count

from gooey import GooeyParser, Gooey

import printers
from distances import K2P_distance, K2P_distance_ambiguity
from functions import read_species, read_fasta, compute_frequencies, remove_duplicates, compute_groups, run_processes, \
    compute_individual_groups


@Gooey(
    program_name="Hooty",
    program_description="Compute distances between groups of sequences",
    use_cmd_args=True,
)
def main():
    start = time.time()

    parser = GooeyParser()
    parser.add_argument(
        "fasta_file",
        metavar="Fasta file",
        help="path to the fasta file",
        widget="FileChooser",
    )
    parser.add_argument(
        "species_file",
        metavar="Species file",
        help="path to the species file",
        widget="FileChooser",
    )
    parser.add_argument(
        "-o",
        "--output",
        metavar="Output file",
        default=None,
        dest="output_file",
    )
    parser.add_argument(
        "-l",
        "--log",
        metavar="Log file",
        default=None,
        dest="log_file",
    )
    parser.add_argument(
        "-p",
        "--processes",
        metavar="Processes",
        default=cpu_count(),
        dest="processes",
        type=int,
        widget="IntegerField",
    )
    parser.add_argument(
        "-u",
        "--unambiguous",
        metavar="Unambiguous",
        default=False,
        action="store_true",
        help="ignore ambiguous bases",
        dest="unambiguous",
    )
    parser.add_argument(
        "-t",
        "--threshold",
        metavar="Threshold",
        default=0.05,
        help="max percentage of ambiguous bases in a group",
        dest="threshold",
        type=float,
    )
    parser.add_argument(
        "--full-matrix",
        metavar="Full matrix",
        default=None,
        help="compute distances for all sequences",
        dest="full_matrix",
    )
    args = parser.parse_args()

    if not os.path.exists(args.fasta_file):
        print("Error: fasta file does not exist")
        return 1
    if not os.path.exists(args.species_file):
        print("Error: species file does not exist")
        return 2

    species, groups, group_offsets = read_species(args.species_file)
    n_groups = len(group_offsets)
    # print(f"Number of species: {len(species)}")
    print(f"Number of groups: {n_groups}")

    dup_seqs = read_fasta(args.fasta_file, species, groups)
    print(f"Number of sequences (with duplicates): {len(dup_seqs)}")

    frequencies = None
    distance_f = K2P_distance
    if not args.unambiguous:
        frequencies = compute_frequencies(dup_seqs, n_groups, args.threshold)
        print("Computed frequencies")
        distance_f = K2P_distance_ambiguity

    seqs, seqs_offsets, min_dist_0 = remove_duplicates(dup_seqs, n_groups)
    print(f"Number of sequences (without duplicates): {len(seqs)}")
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
    return 0


if __name__ == "__main__":
    error_code = main()
    if error_code != 0:
        exit(error_code)
