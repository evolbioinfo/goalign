# Goalign: toolkit and api for alignment manipulation

## Commands

### sw
Aligns 2 sequences using Smith&Waterman algorithm.

Input : Fasta file
Output: Aligned file (format depending on format options)

If neither --match nor --mismatch are specified, then match and mismatch scores
are taken from blosum62 or dnafull substitution matrices (taken from EMBOSS WATER)
depending on the input sequences alphabets.

Score for opening a gap is specified by --gap-open option and score for extending a gap is
specified by --gap-extend option (they should be negative).

Input file must be a fasta file containing 2 sequences. Output format may be specified
by formatting options (-p, -x, etc.).

#### Usage
```
Usage:
  goalign sw [flags]

Flags:
      --gap-extend float   Score for extending a gap  (default -0.5)
      --gap-open float     Score for opening a gap  (default -10)
  -h, --help               help for sw
  -l, --log string         Alignment log file (default "none")
      --match float        Score for a match (if omitted, then take substitution matrix) (default 1)
      --mismatch float     Score for a mismatch (if omitted, then take substitution matrix) (default -1)
  -o, --output string      Alignment output file (default "stdout")

Global Flags:
  -i, --align string    Alignment input file (default "stdin")
      --auto-detect     Auto detects input format (overrides -p, -x and -u)
  -u, --clustal         Alignment is in clustal? default fasta
      --input-strict    Strict phylip input format (only used with -p)
  -x, --nexus           Alignment is in nexus? default fasta
      --no-block        Write Phylip sequences without space separated blocks (only used with -p)
      --one-line        Write Phylip sequences on 1 line (only used with -p)
      --output-strict   Strict phylip output format (only used with -p)
  -p, --phylip          Alignment is in phylip? default fasta
```

#### Examples

seq.fa
```
>nt1
CTGGGGTTTAACCAGCCATGCCAGTGCAGGTTTAAGAACCGATCCGTACTCTGGGTTACTGATGAAGGATGGGCCGTATC
GCCCCCTTGCGACGTTTCCA
>nt2
TATTATCGTATCGTTTGCATAGACCCGTTATGCCAGCAGATACAGCGTCACAAACTTAGGCTGTAGGGCGTTAGCGGCGC
TCCATGTTTAGACTCACGCC
EOF
```

```
goalign sw -i seqs.fa
```

should give:
```
>nt1
GTTT-----AACCAGCCATGCCAGTGCAGGTTTAAGAACCGATCCGT-----ACTCTGGGTTACTGATGAAGGATGGGCC
GTATCGCCCCCTTGCGACGTTTCCA
>nt2
GTTTGCATAGACCCGTTATGCCA--GCAGAT------ACAG---CGTCACAAACTTAGG----CTG--------TAGGGC
GT---------TAGCGGCG-CTCCA
```
