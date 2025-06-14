use std::time::Instant;

use rayon::iter::ParallelIterator;

use crate::structs::Separator;

mod distances;
mod parsing;
mod printers;
mod structs;

use clap::{Parser, arg};
use rayon::prelude::IntoParallelIterator;

use crate::{
	distances::{k2p, k2p_ambiguity},
	parsing::{
		compute_frequencies, get_output_file_path, read_fasta, read_species, remove_duplicates,
	},
};

#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
struct Args {
	#[arg(help = "path to the fasta file", name = "fasta file")]
	fasta_file: String,
	#[arg(help = "path to the species file", name = "species file")]
	species_file: String,
	#[arg(
		help = "path to the output file",
		name = "output file",
		long = "output",
		short
	)]
	output_file: Option<String>,
	#[arg(
		help = "path to the full matrix file",
		name = "full matrix file",
		long = "full_matrix",
		short
	)]
	full_matrix_file: Option<String>,
	#[arg(
		help = "separator to use in the output file",
		long,
		short,
		default_value = "semicolon"
	)]
	separator: Separator,
	#[arg(
		help = "max percentage of ambiguous bases in a group",
		long,
		short,
		default_value_t = 0.0
	)]
	threshold: f32,
	#[arg(help = "treat ambiguous sites as similarities", long, short)]
	unambiguous: bool,
}

fn main() {
	let start = Instant::now();
	let args = Args::parse();

	let (species, groups, group_offsets) = read_species(args.species_file);
	let n_groups = group_offsets.len();
	println!("Number of groups: {n_groups}");

	let dup_seqs = read_fasta(&args.fasta_file, &species, &groups);
	println!("Number of sequences (with duplicates): {}", dup_seqs.len());

	let (frequencies, distance_f): (_, fn(_, _, _) -> _) = if args.unambiguous {
		(None, k2p)
	} else {
		let frequencies = Some(compute_frequencies(&dup_seqs, n_groups, args.threshold));
		println!("Computed frequencies");
		(frequencies, k2p_ambiguity)
	};

	let (seqs, seqs_offsets, min_dist_0) = remove_duplicates(&dup_seqs, n_groups);
	println!("Number of sequences (without duplicates): {}", seqs.len());

	let mut result: Vec<_> = (0..n_groups)
		.map(|i| vec![(f64::NAN, f64::NAN); i + 1])
		.collect();
	let transposed: Vec<_> = (0..n_groups)
		.into_par_iter()
		.map(|g1| {
			let mut result = Vec::with_capacity(n_groups - g1);
			let g1_offset = seqs_offsets[g1];
			#[allow(clippy::needless_range_loop)]
			for g2 in g1..n_groups {
				let g2_offset = seqs_offsets[g2];
				let (mut min, mut max) = (f64::INFINITY, f64::NEG_INFINITY);
				for i in g1_offset.offset..g1_offset.offset + g1_offset.count {
					for j in g2_offset.offset..g2_offset.offset + g2_offset.count {
						if i == j {
							continue;
						}
						let distance = distance_f(&seqs[i], &seqs[j], frequencies.as_ref());
						if distance < min {
							min = distance;
						}
						if distance > max {
							max = distance;
						}
					}
				}
				if g1 == g2 {
					if min_dist_0[g1] {
						min = 0.0;
						if g1_offset.count == 1 {
							max = 0.0;
						}
					} else if g1_offset.count == 1 {
						max = f64::INFINITY;
						min = f64::NEG_INFINITY;
					}
				}
				result.push((min, max));
			}
			result
		})
		.collect();

	for g1 in 0..transposed.len() {
		for g2 in 0..transposed[g1].len() {
			result[g2 + g1][g1] = transposed[g1][g2];
		}
	}
	let result_str = printers::to_sv(&result, &species, &groups, args.separator, 2);
	let dest_path = args
		.output_file
		.unwrap_or(get_output_file_path(&args.fasta_file));
	if let Err(e) = std::fs::write(&dest_path, result_str) {
		eprintln!("An error occurred while writing to output: {e}")
	} else {
		println!("Output written to {dest_path}");
	}

	if let Some(full_matrix_path) = args.full_matrix_file {
		let n_groups = dup_seqs.len();
		let mut full_result: Vec<_> = (0..n_groups).map(|i| vec![f64::NAN; i + 1]).collect();
		let transposed: Vec<_> = (0..n_groups)
			.into_par_iter()
			.map(|g1| {
				let mut result = Vec::with_capacity(n_groups - g1);
				for g2 in g1..n_groups {
					result.push(distance_f(
						&dup_seqs[g1],
						&dup_seqs[g2],
						frequencies.as_ref(),
					));
				}
				result
			})
			.collect();

		for g1 in 0..transposed.len() {
			for g2 in 0..transposed[g1].len() {
				full_result[g2 + g1][g1] = transposed[g1][g2];
			}
		}
		let full_result_str = printers::to_sv(
			&full_result,
			&dup_seqs.iter().map(|s| s.name.clone()).collect::<Vec<_>>(),
			&(0..n_groups).collect(),
			args.separator,
			15,
		);
		if let Err(e) = std::fs::write(&full_matrix_path, full_result_str) {
			eprintln!("An error occurred while writing full matrix: {e}")
		} else {
			println!("Full matrix written to {full_matrix_path}");
		}
	}
	println!("Completion time: {:?}", Instant::now() - start)
}
