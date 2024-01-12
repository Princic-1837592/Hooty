# Examples

In this directory you can find a collection of example input and output files to use with Hooty.\
In each example you will find:

- `alignment.fasta`: the file containing DNA sequences
- `partition.txt`: the file containing the partitioning of the sequences into groups
- `expected.csv`: the expected output file
- `expected_fullmatrix.csv`: the expected output file with the full matrix

To run Hooty on an example, you can use the following command:

```bash
python hooty.py examples/example1/alignment.fasta examples/example1/partition.txt
```

This will create the output file `alignment.csv` (maintaining the name of the fasta file) in the same directory as the
input files.\
Alternatively, you can specify the output file name and path with the `-o` option:

```bash
python hooty.py examples/example1/alignment.fasta examples/example1/partition.txt -o examples/example1/output.csv
```

You can also run the example with the option `--full-matrix` to get the full pairwise matrix, specifying the output file:

```bash
python hooty.py examples/example1/alignment.fasta examples/example1/partition.txt --full-matrix examples/example1/fullmatrix.csv
```

More information about the options can be found on
the [documentation](https://princic-1837592.github.io/Hooty/index.html#options)
