use crate::structs::Separator;

#[allow(private_bounds)]
pub(crate) fn to_sv(
	matrix: &Vec<Vec<impl Cell>>,
	species: &Vec<String>,
	groups: &Vec<usize>,
	separator: Separator,
	precision: usize,
) -> String {
	let separator = separator.to_string();
	let mut header_groups = vec![vec![]; groups.last().unwrap() + 2];
	for (species, group) in species.iter().zip(groups) {
		header_groups[group + 1].push(species.clone());
	}
	let header: Vec<_> = header_groups.iter().map(|hg| hg.join("-")).collect();
	let mut lines = Vec::with_capacity(1 + matrix.len());
	lines.push(header.join(separator));
	for (row, group) in matrix.iter().zip(header.iter().skip(1)) {
		let numbers = row
			.iter()
			.map(|c| c.format(precision))
			.collect::<Vec<_>>()
			.join(separator);
		lines.push(format!("{}{}{}", group, separator, numbers));
	}
	lines.join("\n")
}

trait Cell {
	fn format(&self, precision: usize) -> String;
}

impl Cell for (f64, f64) {
	fn format(&self, precision: usize) -> String {
		format!(
			"{} - {}",
			self.0.format(precision),
			self.1.format(precision)
		)
	}
}

impl Cell for f64 {
	fn format(&self, precision: usize) -> String {
		if self.abs() == 0.0 {
			"0".to_owned()
		} else if self.is_nan() || self.is_infinite() {
			"/".to_owned()
		} else {
			format!("{:.precision$}", self * 100.0)
		}
	}
}
