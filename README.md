# Hooty

Hooty computes the minimum and maximum K2P distance between predetermined groups
in an alignment and outputs results in a matrix format.
It also includes an option for the interpretation of ambiguous bases.

## Usage

Hooty takes as input a FASTA file and a text file containing a partition of species.
The output is a matrix with the minimum and maximum K2P distances between the groups.

```
cargo install hooty
```

```
hooty <alignment.fasta> <partition.txt>
```

Hooty can also compute the full pairwise K2P distances between all sequences in the alignment.

```
hooty <alignment.fasta> <partition.txt> --full-matrix <full_matrix.csv>
```

## Documentation

You can find the complete documentation and all the available options for Hooty [here](https://princic-1837592.github.io/Hooty/).
