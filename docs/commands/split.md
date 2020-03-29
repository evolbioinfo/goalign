# Goalign: toolkit and api for alignment manipulation

## Commands

### split
This command splits an input alignment according to partitions given as input.

The partitions are defined as in [RAxML](https://cme.h-its.org/exelixis/web/software/raxml/index.html).

#### Usage
```
goalign split -i align.phylip --partition partition.txt


Usage:
  goalign split [flags]

Flags:
  -h, --help                help for split
  -o, --out-prefix string   Prefix of output files
      --partition string    File containing definition of the partitions (default "none")

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

* Spliting an alignment

input.fa
```
>s1
AAAACCCCCGG
>2
AAAACCCCCGG
>3
AAAACCCCCGG
>4
AAAACCCCCGG
>5
AAAACCCCCGG
```

partition.txt
```
M1,p1=1-4,10-11
M2,p2=5-9
```

This command:
```
goalign split -i input.fa --partition partition.txt --out-prefix ./
```

Should produce:

p1.fa
```
>s1
AAAAGG
>2
AAAAGG
>3
AAAAGG
>4
AAAAGG
>5
AAAAGG
```

and p2.fa
```
>s1
CCCCC
>2
CCCCC
>3
CCCCC
>4
CCCCC
>5
CCCCC
```
