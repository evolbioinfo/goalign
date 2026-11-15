# Goalign: toolkit and api for alignment manipulation

## Commands

### compute
This command implements different computations:
1. `goalign compute distance` : Computes a distance matrix from an input DNA alignment, with different evolutionary models. In the case of ambigous nucleotides (IUPAC), one mutation is counted if characters are not compatible (ex: R vs. Y). Possible models are:
    - pdist
	  - rawdist : Raw distance (like pdist, without normalization by length)
    - jc      : Juke-Cantor
    - k2p     : Kimura 2 Parameters
    - f81     : Felsenstein 81
    - f84     : Felsenstein 84
    - tn93    : Tamura and Nei 1993
  If distance is pdist (nucleotides), then giving the option --rm-ambiguous will not take into 
  account ambiguous positions that compatible, for length normalization.
  For example if --rm-ambiguous is given, then R vs. Y will be taken into account
  because there is a difference. And N vs. A won't be taken into account in total length
  because we are not sure whether they are identical.
  In case of a nucleotidic alignment, it is possible to specify sequence ranges to compare. For example, 
  goalign compute distance -m pdist -i align.ph -p --range1 0:9 --range2 10:19
  will compute distance only between sequences [0 to 9] and sequences [10 to 19].
  Output matrix will be formatted the same way as usual, except that it will be made of 0 except for
  the comparisons 0 vs. 10; 0 .vs 11; ...; 9 vs. 19.
2. `goalign compute entropy`: Computes the Shannon entropy of each sites of the input alignment or the average entropy of all sites (`-a` option).
    - Entropy is computed with the natural logarithm by default, i.e. values are in **nats**. `--log-base 2` gives values in bits, `--log-base 10` in base 10 (any positive base is accepted).
    - Characters are upper-cased before counting. `*` and `.` are never taken into account; `-` is not taken into account with `--remove-gaps` (`-g`); N, X, IUPAC ambiguity codes and any other non standard nucleotide/amino acid character are not taken into account with `--ignore-ambiguities`. Without this option, they are counted as states of their own, which may strongly inflate entropy on alignments with many undetermined positions (e.g. viral genomes).
    - With `--ignore-ambiguities`, the alphabet matters: an alignment auto-detected as "both" (e.g. a protein alignment made only of A, C, G, T, N...) is treated as nucleotides (use `--alphabet aa` to force amino acids), and an alignment of unknown alphabet (characters that are neither nucleotides nor amino acids) produces an error.
    - Sites with no character taken into account have an entropy of `NaN`, and are not taken into account in the average.
    - Values are printed with 3 decimals by default; use `--precision` to change it (applies to both per-site values and the `-a` average). Low entropy values (e.g. a single variant among thousands of sequences) need more than 3 decimals to be meaningful.
3. `goalign compute pssm`: Computes and prints a Position specific scoring matrix. Different kind of matrices may be computed, depending on `-n` option:
    - `-n 0` : None, means raw counts
    - `-n 1` : By column frequency, i.e. frequency of nt/aa per site/column
    - `-n 2` : By column frequency compared to alignment frequency: same as -n 1, but divides by frequency of the nt/aa in the whole alignment
    - `-n 3` : By column frequency compared to uniform frequency: same as -n 1, but divides by uniform frequency of the nt/aa (1/4 for nt, 1/20 for aa)
    - `-n 4` : Normalization "Logo". Uses Shannon entropy in **bits** (log2), as is the convention for sequence logos (unlike `goalign compute entropy`, which outputs nats by default).
	Option `-c` allows to add pseudo counts before normalization, and option `-l` log2 transforms the values.
 4. `goalign compute simplot`: See the [dedicated page](simplot.md). Compute a similarity plot between a query sequence and other sequences, using a sliding window.

#### Usage

* General command
```
Usage:
  goalign compute [command]

Available Commands:
  distance    Compute distance matrix from an input alignment
  entropy     Computes entropy of a given alignment
  pssm        Computes and prints a Position specific scoring matrix
  [simplot](simplot.md))     Computes simplot data and image

Flags:
  -h, --help   help for compute

Global Flags:
  -i, --align string   Alignment input file (default "stdin")
  -p, --phylip         Alignment is in phylip? False=Fasta
  --input-strict       Strict phylip input format (only used with -p)
```

* entropy command
```
goalign compute entropy [flags]

Flags:
  -a, --average              Compute only the average entropy of input alignment
  -h, --help                 help for entropy
      --ignore-ambiguities   If true, then do not take into account N, X and other ambiguous/non standard characters in the computation
      --log-base string      Base of the logarithm: "e" (nats, default), 2 (bits), 10, or any positive number (default "e")
      --precision int        Number of decimals of output entropy values (per site and average) (default 3)
  -g, --remove-gaps          If true, then do not take into account gaps in the computation
```

* distance command
```
goalign compute distance [flags]

Flags:
  -a, --average         Compute only the average distance between all pairs of sequences
  -m, --model string    Model for distance computation (default "k2p")
  -o, --output string   Distance matrix output file (default "stdout")
  -r, --rm-gaps         Do not take into account positions containing >=1 gaps

Global Flags:
  -i, --align string   Alignment input file (default "stdin")
  -p, --phylip         Alignment is in phylip? False=Fasta
  --input-strict       Strict phylip input format (only used with -p)
```

* pssm command
```
Usage:
  goalign compute pssm [flags]

Flags:
  -l, --log                   (normalized) Values in log2
  -n, --normalization int     Counts normalization
  -c, --pseudo-counts float   Value added to (normalized) counts

Global Flags:
  -i, --align string   Alignment input file (default "stdin")
  -p, --phylip         Alignment is in phylip? False=Fasta
```

#### Examples

* Generating a random tree with 5 tips ([Gotree](https://github.com/evolbioinfo/gotree)), simulating an alignment from this tree ([seq-gen](https://github.com/rambaut/Seq-Gen), and computing a distance matrix (model f81) from this alignment:
```
gotree generate yuletree -l 5 --seed 1 -o true_tree.nw
seq-gen -op -mGTR -l500 -z 2 -n 1 true_tree.nw > alignment.phy
goalign compute distance -i alignment.phy -m f81 -o dist.txt -p -t 10
```

Should give the following distance matrix:

```
5
Tip4    0.000000000000  0.174911845895  0.192803956978  0.232646053483  0.235379041630
Tip0    0.174911845895  0.000000000000  0.082364641962  0.128396525775  0.142776083476
Tip3    0.192803956978  0.082364641962  0.000000000000  0.071285264523  0.086842665158
Tip2    0.232646053483  0.128396525775  0.071285264523  0.000000000000  0.111961817720
Tip1    0.235379041630  0.142776083476  0.086842665158  0.111961817720  0.000000000000
```

* Generating a random tree with 100 tips ([Gotree](https://github.com/evolbioinfo/gotree)), simulating an alignment from this tree ([seq-gen](https://github.com/rambaut/Seq-Gen), and computing entropy of each site of this alignment:
```
gotree generate yuletree -l 200 --seed 1 -o true_tree.nw
seq-gen -op -mGTR -l10 -z 2 -n 1 true_tree.nw > alignment.phy
goalign compute entropy  -i alignment.phy -p 
```

Should give the following matrix:
```
Alignment  Site  Entropy
0          0     0.616
0          1     0.858
0          2     0.806
0          3     0.621
0          4     0.621
0          5     0.758
0          6     0.813
0          7     0.780
0          8     0.762
0          9     0.816
```

* Generating a random tree with 100 tips ([Gotree](https://github.com/evolbioinfo/gotree)), simulating an alignment from this tree ([seq-gen](https://github.com/rambaut/Seq-Gen)), and computing a logo from this alignment:
```
gotree generate yuletree -l 200 --seed 1 -o true_tree.nw
seq-gen -op -mGTR -l10 -z 2 -n 1 true_tree.nw > alignment.phy
goalign compute pssm -n 4 -i alignment.phy -p -c 0.0001
```

Should give the following matrix:
```
    A      C      G      T
1   0.089  0.901  0.122  0.000
2   0.541  0.038  0.038  0.145
3   0.628  0.033  0.117  0.059
4   0.088  0.022  0.077  0.916
5   0.022  0.077  0.916  0.088
6   0.118  0.063  0.027  0.698
7   0.066  0.628  0.058  0.074
8   0.096  0.674  0.061  0.044
9   0.054  0.703  0.054  0.090
10  0.041  0.576  0.189  0.016
```
