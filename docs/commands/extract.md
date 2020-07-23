# Goalign: toolkit and api for alignment manipulation

## Commands

### extract
This command extracts several sub-alignments from an input alignment. It is similar to [subseq](subseq.md), with two main differences:

1. As input, it takes an annotation file defining the coordinates of sub-alignments to extract, and can then extract several sub-alignments in one command;
2. Each sub-alignment may be defined as several "blocks" (~bed format blocks), potentially overlapping (even in any order). This may be useful if sequences are subject to ~frameshifts (see [this case](https://www.ncbi.nlm.nih.gov/protein/YP_009724389.1?from=4393&to=5324) for example).

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
  goalign extract [flags]

Flags:
      --coordinates string   File with all coordinates of the sequences to extract (default "none")
  -h, --help                 help for extract
  -o, --output string        Output folder (default ".")
      --ref-seq string       Reference sequence on which coordinates are given (default "none")
      --translate int        Wether the extracted sequence will be translated (only if input alignment is nucleotide). <0: No translation, 0: Std code, 1: Vertebrate mito, 2: Invertebrate mito (default -1)

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

* Extract 2 sub sequences from an input alignment:

Input alignment `al.fa`:
```
>s1
ACGTACGT
>s2
-CGT-C-T
>s3
ACGTACGT
>s4
ACGTTCGA
>s5
ACGTTCGA
```

Annotation file `coords.txt`:
```
1,2	4,5	output.1
3	6	output.2
```

This annotation file defines two sub-alignments:

1. A sub-alignment with two overlapping blocks: `[1,4[+[2,5[` (or `[1,3]+[2,4]`)
2. A sub-alignment with 1 block: `[3,6[` (or `[3,5]`)


Command:
```
> goalign extract -i al.fa --coordinates coords.txt -o out
```

Should produce 2 files in the `out` directory:

`out/output.1.fa`:
```
>s1
CGTGTA
>s2
CGTGT-
>s3
CGTGTA
>s4
CGTGTT
>s5
CGTGTT
```

`out/output.2.fa`:
```
>s1
TAC
>s2
T-C
>s3
TAC
>s4
TTC
>s5
TTC
EOF
```
