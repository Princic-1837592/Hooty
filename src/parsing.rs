use crate::structs::{AmbiguityInfo, Base, Frequencies, Offset, Sequence, REPLACEMENTS};
use std::collections::HashSet;
use std::path::Path;

pub(crate) fn read_species<P: AsRef<Path>>(path: P) -> (Vec<String>, Vec<usize>, Vec<usize>) {
	let mut species = Vec::new();
	let mut groups = Vec::new();
	let mut offsets = Vec::new();
	let file = std::fs::read_to_string(path).expect("Failed to read species file");
	let mut group = 0;
	for line in file.lines() {
		offsets.push(species.len());
		for s in line.split(',').map(&str::trim) {
			groups.push(group);
			species.push(s.to_owned());
		}
		group += 1;
	}
	(species, groups, offsets)
}

fn to_single_lines<'a>(mut file: String) -> Vec<String> {
	file.push_str("\n>");
	let mut buffer = "".to_owned();
	let mut result = Vec::from([file
		.lines()
		.next()
		.expect("Expected line. If you see this please contact the developer")
		.to_owned()]);
	for line in file.lines().skip(1) {
		if line.starts_with('>') {
			result.push(buffer.clone());
			result.push(line.to_owned());
			buffer.clear()
		} else {
			buffer.push_str(line);
		}
	}
	result.pop();
	result
}

pub(crate) fn read_fasta<P: AsRef<Path>>(
	path: P,
	species: &Vec<String>,
	groups: &Vec<usize>,
) -> Vec<Sequence> {
	let file = std::fs::read_to_string(path).expect("Failed to read fasta file");
	let first_three: Vec<_> = file.lines().filter(|l| l.len() > 0).take(3).collect();
	let lines: Vec<_> = if first_three.len() > 2
		&& first_three[2]
			.chars()
			.next()
			.expect("Line should not be empty. If you see this please contact the developer")
			!= '>'
	{
		to_single_lines(file)
	} else {
		file.lines()
			.filter_map(|l| (l.len() > 0).then_some(l.to_owned()))
			.collect()
	};

	let mut sequences = Vec::new();
	for l in (0..lines.len()).step_by(2) {
		let name = lines[l][1..].trim().to_owned();
		let of_groups = species
			.iter()
			.enumerate()
			.filter_map(|(i, s)| name.contains(s).then_some(groups[i]))
			.collect::<HashSet<_>>();
		match of_groups.len() {
			1 => sequences.push(Sequence {
				index: l / 2,
				group: *of_groups.iter().next().unwrap(),
				name,
				seq: lines[l + 1].chars().map(Base::from).collect(),
			}),
			0 => {
				eprintln!("WARNING: match not found for sequence '{}'", name)
			}
			_ => {
				eprintln!(
					"WARNING: species of sequence '{}' is not unique. Found matches for {} species",
					name,
					of_groups.len()
				)
			}
		}
	}
	sequences
}

pub(crate) fn compute_frequencies(
	sequences: &Vec<Sequence>,
	n_groups: usize,
	threshold: f32,
) -> AmbiguityInfo {
	let mut groups = vec![vec![Frequencies::new(); sequences[0].seq.len()]; n_groups];
	for seq in sequences {
		for (i, &base) in seq.seq.iter().enumerate() {
			if base <= Base::T {
				groups[seq.group][i].count[base as usize] += 1;
				groups[seq.group][i].normal_count += 1;
			} else {
				groups[seq.group][i].ambiguous_count += 1
			}
		}
	}
	for group in &mut groups {
		for freq in group {
			for (i, replace) in REPLACEMENTS.iter().copied().enumerate() {
				let valid_replace = replace.iter().filter(|&&b| (b != Base::NoneBase));
				let rep_count = 1.max(valid_replace.clone().map(|&b| freq.count[b as usize]).sum());
				for b in valid_replace.copied() {
					freq.frequencies[i][b as usize] =
						freq.count[b as usize] as f64 / rep_count as f64;
				}
			}
			freq.percentage = freq.ambiguous_count as f32
				/ 1.max(freq.ambiguous_count + freq.normal_count) as f32;
		}
	}
	AmbiguityInfo { groups, threshold }
}

pub(crate) fn remove_duplicates(
	sequences: &Vec<Sequence>,
	n_groups: usize,
) -> (Vec<Sequence>, Vec<Offset>, Vec<bool>) {
	let mut min_dist_0 = vec![false; n_groups];
	let mut dedup = HashSet::new();
	for seq in sequences {
		if !dedup.contains(seq) {
			dedup.insert(seq.clone());
		} else {
			min_dist_0[seq.group] = true;
		}
	}
	let mut sequences: Vec<_> = dedup.into_iter().collect();
	sequences.sort_by_key(|s| (s.group, s.index));
	let mut offsets = vec![
		Offset {
			offset: usize::MAX,
			count: 0
		};
		n_groups
	];
	for (i, seq) in sequences.iter().enumerate() {
		if offsets[seq.group].offset == usize::MAX {
			offsets[seq.group].offset = i;
		}
		offsets[seq.group].count += 1;
	}
	(sequences, offsets, min_dist_0)
}

pub(crate) fn get_output_file_path(path: &String) -> String {
	let path = Path::new(path);
	let stem = path.file_stem().and_then(|s| s.to_str()).unwrap_or("");
	let parent = path.parent().unwrap_or(Path::new(""));

	format!(
		"{}.csv",
		if parent.as_os_str().is_empty() {
			stem.to_string()
		} else {
			parent.join(stem).to_string_lossy().into_owned()
		}
	)
}
