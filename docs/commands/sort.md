# Goalign: toolkit and api for alignment manipulation

## Commands

### sort
This command sort an input alignment by sequence names.

#### Usage
```
Usage:
  goalign sort [flags]

Flags:
  -o, --output string   Sorted alignment output file (default "stdout")

Global Flags:
  -i, --align string    Alignment input file (default "stdin")
  --input-strict    Strict phylip input format (only used with -p)
  -x, --nexus           Alignment is in nexus? default fasta
  --output-strict   Strict phylip output format (only used with -p)
  -p, --phylip          Alignment is in phylip? default fasta
```

#### Examples

* Sorting a random alignment
```
goalign random -s 10 -l 10 -n 5 | goalign shuffle seqs | goalign sort
```

It should give the following alignment:
```
>Seq0000
GATTAATTTG
>Seq0001
CCGTAGGCCA
>Seq0002
GAATCTGAAG
>Seq0003
ATCGAACACT
>Seq0004
TTAAGTTTTC
```
