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

import os
import time
from math import inf

import printers
from distances import K2P_distance, K2P_distance_ambiguity
from functions import read_species, read_fasta, compute_frequencies, remove_duplicates, compute_groups, run_processes, \
    compute_individual_groups


def main(parser) -> int:
    args = parser.parse_args()
    start = time.time()

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

    string = printers.to_sv(result, species, groups, args.separator, 2)
    if args.output_file is None:
        args.output_file = os.path.splitext(args.fasta_file)[0] + ".csv"
    with open(args.output_file, "w", encoding="utf-8") as f:
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
        string = printers.to_sv(
            full_result,
            [seq.name for seq in dup_seqs],
            [i for i in range(n_groups)],
            args.separator,
            15
        )
        with open(args.full_matrix, "w", encoding="utf-8") as f:
            f.write(string)
        print(f"Full matrix written to {args.full_matrix}")

    end = time.time()
    print(f"{end - start:.2f}")
    return 0
