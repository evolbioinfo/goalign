# Goalign: toolkit and api for alignment manipulation

## Commands

### identical
Tells wether two alignments are identical.

Alignments are considered identical if:

1. They have the same number of sequences;
2. Each sequence of the first alignment have a corresponding sequence
   in the second alignment having the same name and the same sequence.

Identical alignments may have sequences in different order

Example:


align1.fa

```
>Seq0000
GATTAA---GCCGTAGGCCAGAATCTGAAG
>Seq0001
ATCGAA---TTTAAGTTT---CTTCTAATG
>Seq0002
GAGAGGACTAGTTCATACTTTTTAAACACT
```

align2.fa

```
>Seq0001
ATCGAA---TTTAAGTTT---CTTCTAATG
>Seq0000
GATTAA---GCCGTAGGCCAGAATCTGAAG
>Seq0002
GAGAGGACTAGTTCATACTTTTTAAACACT
```

```
goalign identical -i align1.fa -c align2.fa
```

should print:

```
true
```


#### Usage
```
Usage:
  goalign identical [flags]

Flags:
  -c, --compared string   Compared alignment file (default "none")

Global Flags:
  -i, --align string    Alignment input file (default "stdin")
      --auto-detect     Auto detects input format (overrides -p, -x and -u)
  -u, --clustal         Alignment is in clustal? default fasta
      --input-strict    Strict phylip input format (only used with -p)
  -x, --nexus           Alignment is in nexus? default fasta
      --output-strict   Strict phylip output format (only used with -p)
  -p, --phylip          Alignment is in phylip? default fasta
```
