use clap::ValueEnum;
use std::hash::{Hash, Hasher};
use std::ops::Index;
use Base::*;

#[derive(Eq, Debug, Clone)]
pub(crate) struct Sequence {
	pub(crate) index: usize,
	pub(crate) group: usize,
	pub(crate) name: String,
	pub(crate) seq: Vec<Base>,
}

impl Hash for Sequence {
	fn hash<H: Hasher>(&self, state: &mut H) {
		self.group.hash(state);
		self.seq.hash(state);
	}
}

impl PartialEq for Sequence {
	fn eq(&self, other: &Self) -> bool {
		self.group.eq(&other.group) && self.seq.eq(&other.seq)
	}
}

#[derive(Copy, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug)]
pub(crate) enum Base {
	A,
	C,
	G,
	T,
	W,
	S,
	M,
	K,
	R,
	Y,
	B,
	D,
	H,
	V,
	N,
	GAP,
	NoneBase,
}

impl From<char> for Base {
	fn from(value: char) -> Self {
		match value {
			'A' => A,
			'C' => C,
			'G' => G,
			'T' => T,
			'W' => W,
			'S' => S,
			'M' => M,
			'K' => K,
			'R' => R,
			'Y' => Y,
			'B' => B,
			'D' => D,
			'H' => H,
			'V' => V,
			'N' => N,
			_gap => GAP,
		}
	}
}

pub(crate) const REPLACEMENTS: [[Base; 3]; 10] = [
	[A, T, NoneBase],
	[G, C, NoneBase],
	[A, C, NoneBase],
	[G, T, NoneBase],
	[A, G, NoneBase],
	[C, T, NoneBase],
	[C, G, T],
	[A, G, T],
	[A, C, T],
	[A, C, G],
];

#[derive(Debug, Clone, ValueEnum, Copy)]
pub(crate) enum Separator {
	Comma,
	Semicolon,
	Tab,
}

impl Separator {
	pub(crate) fn to_string<'a>(self) -> &'a str {
		match self {
			Separator::Comma => ",",
			Separator::Semicolon => ";",
			Separator::Tab => "\t",
		}
	}
}

#[derive(Clone)]
pub(crate) struct Frequencies {
	pub(crate) count: [usize; 4],
	pub(crate) frequencies: [[f64; 4]; 10], //todo REPLACEMENTS
	pub(crate) normal_count: usize,
	pub(crate) ambiguous_count: usize,
	pub(crate) percentage: f32,
}

impl Frequencies {
	pub(crate) fn new() -> Self {
		Self {
			count: [0, 0, 0, 0],
			frequencies: [[0.0; 4]; 10],
			normal_count: 0,
			ambiguous_count: 0,
			percentage: 0.0,
		}
	}
}

impl Index<usize> for Frequencies {
	type Output = [f64; 4];

	fn index(&self, index: usize) -> &Self::Output {
		&self.frequencies[index]
	}
}

pub(crate) struct AmbiguityInfo {
	pub(crate) groups: Vec<Vec<Frequencies>>,
	pub(crate) threshold: f32,
}

#[derive(Copy, Clone, Debug)]
pub(crate) struct Offset {
	pub(crate) offset: usize,
	pub(crate) count: usize,
}
