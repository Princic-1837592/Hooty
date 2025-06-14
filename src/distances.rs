use crate::structs::{AmbiguityInfo, Base::*, Sequence};

const TRANSITIONS: [[u8; 4]; 4] = [[0, 0, 1, 0], [0, 0, 0, 1], [1, 0, 0, 0], [0, 1, 0, 0]];

pub(crate) fn k2p(seq1: &Sequence, seq2: &Sequence, _: Option<&AmbiguityInfo>) -> f64 {
	let pairs = seq1.seq.iter().zip(&seq2.seq).filter(|&(&b1, &b2)| {
		b1 != Gap && b1 != N && b1 != NoneBase && b2 != Gap && b2 != N && b2 != NoneBase
	});
	let (mut ts_count, mut tv_count, length) = (0, 0, pairs.clone().count());
	if length == 0 {
		return f64::NAN;
	}
	for (&b1, &b2) in pairs.filter(|&(&b1, &b2)| b1 != b2 && b1 <= T && b2 <= T) {
		ts_count += TRANSITIONS[b1 as usize][b2 as usize];
		tv_count += 1 - TRANSITIONS[b1 as usize][b2 as usize];
	}
	let p = ts_count as f64 / length.max(1) as f64;
	let q = tv_count as f64 / length.max(1) as f64;
	-0.5 * ((1.0 - 2.0 * p - q) * (1.0 - 2.0 * q).sqrt()).ln()
}

pub(crate) fn k2p_ambiguity(
	seq1: &Sequence,
	seq2: &Sequence,
	frequencies: Option<&AmbiguityInfo>,
) -> f64 {
	let frequencies = frequencies
		.expect("Expected frequencies in k2p_ambiguity. If you see this contact the developer");
	let pairs = seq1
		.seq
		.iter()
		.zip(&seq2.seq)
		.enumerate()
		.filter(|&(_, (&b1, &b2))| {
			b1 != Gap && b1 != N && b1 != NoneBase && b2 != Gap && b2 != N && b2 != NoneBase
		});
	let (mut ts_count, mut ts_count_f, mut tv_count, mut tv_count_f, mut length) =
		(0, 0.0, 0, 0.0, pairs.clone().count());
	let (g1, g2) = (
		&frequencies.groups[seq1.group],
		&frequencies.groups[seq2.group],
	);

	for (i, (&b1, &b2)) in pairs.filter(|(_, (b1, b2))| b1 != b2) {
		if b1 <= T && b2 <= T {
			ts_count += TRANSITIONS[b1 as usize][b2 as usize];
			tv_count += 1 - TRANSITIONS[b1 as usize][b2 as usize];
		} else if g1[i].percentage < frequencies.threshold
			&& g2[i].percentage < frequencies.threshold
		{
			let b1 = g1[i][b1 as usize - W as usize];
			let b2 = g2[i][b2 as usize - W as usize];
			ts_count_f += b1[A as usize] * b2[G as usize]
				+ b1[G as usize] * b2[A as usize]
				+ b1[C as usize] * b2[T as usize]
				+ b1[T as usize] * b2[C as usize];
			tv_count_f += b1[A as usize] * b2[T as usize]
				+ b1[T as usize] * b2[A as usize]
				+ b1[A as usize] * b2[C as usize]
				+ b1[C as usize] * b2[A as usize]
				+ b1[G as usize] * b2[T as usize]
				+ b1[T as usize] * b2[G as usize]
				+ b1[G as usize] * b2[C as usize]
				+ b1[C as usize] * b2[G as usize];
		} else {
			length -= 1;
		}
	}
	if length == 0 {
		f64::NAN
	} else {
		let p = (ts_count as f64 + ts_count_f) / length as f64;
		let q = (tv_count as f64 + tv_count_f) / length as f64;
		-0.5 * ((1.0 - 2.0 * p - q) * (1.0 - 2.0 * q).sqrt()).ln()
	}
}
