# Goalign: toolkit and api for alignment manipulation

## Commands

### extract
This command extracts several sub-alignments from an input alignment. It is similar to [subseq](subseq.md), with two main differences:

1. As input, it takes an annotation file defining the coordinates of sub-alignments to extract, and can then extract several sub-alignments in one command;
2. Each sub-alignment may be defined as several "blocks" (~bed format blocks), potentially overlapping (even in any order).

`extract` takes an alignment and extracts several sub-alignments from it. Subs-alignments are defined in an input tab separated file with the following mandatory columns:

1. Start coordinates: 0-based inclusive. If the sub-alignment is defined by several blocks, several start coordinates may be given, and separated by comas;
2. End coordinates: 1-based (or 0-based exclusive). If the sub-alignment is defined by several blocks, several end coordinates may be given (same number as start coordinates), and coma separated;
3. Name of the subsequence

Example of an annotation file:

```
0	10	orf1
10	100	orf2
100,105	106,110	orf3
```

The 3rd line defines a sub-alignment containing positions `[100-106[+[105-110[` (or `[100-105]+[105-109]`).

If start (or end) coordinates are outside the alignment, or are not compatible (`start>=end`)then it exits with an error.

If a sub-alignment is defined by several blocks, they are allowed to overlap or be in any order.

Output file names will be defined by the names of the subsequences, with .fa or .phy extension depending on the input file format.

If --ref-seq is given, then the coordinates are defined wrt the given reference sequence (gaps are not taken into acount although they are still present in the output sub-alignment).

If --translate is >=0 and the input alignment is nucleotidic, extracted subsequences are translated into amino acids.
- If --translate < 0 : No translation
- If --translate 0: Standard genetic code
- If --translate 1: Vertebrate mitochondrial genetic code
- If --translate 2: Invertebrate mitochondrial genetic code


For example:
goalign extract -i alignment.fasta -f annotations.txt

If the input file contains several alignments, only the first one is considered.

#### Usage
```
Usage:
  goalign consensus [flags]

Flags:
  --exclude-gaps        Exclude gaps in the majority computation
  -h, --help            help for consensus
  -o, --output string   Alignment output file (default "stdout")

Global Flags:
  -i, --align string       Alignment input file (default "stdin")
      --auto-detect        Auto detects input format (overrides -p, -x and -u)
  -u, --clustal            Alignment is in clustal? default fasta
      --ignore-identical   Ignore duplicated sequences that have the same name and same sequences
      --input-strict       Strict phylip input format (only used with -p)
  -x, --nexus              Alignment is in nexus? default fasta
      --no-block           Write Phylip sequences without space separated blocks (only used with -p)
      --one-line           Write Phylip sequences on 1 line (only used with -p)
      --output-strict      Strict phylip output format (only used with -p)
  -p, --phylip             Alignment is in phylip? default fasta
```

#### Examples

* Consensus of 3 sequences:

Input alignment `al.fa`:
```
>s1
ACGACGACGACC
>2
ATCTT-TTTTTC
>3
ATCTT-TTTTTT
```

```
$ goalign consensus -i al.fa

>consensus
ATCTT-TTTTTC
```
