# Goalign: toolkit and api for alignment manipulation

## Commands

### clean
This command removes alignment sites constitued of >= than a given proportion of gaps. Exception for a cutoff of 0: removes sites constitued of > 0 gaps.
Examples:
- With a cutoff of 0.5: a site with 5 gaps over 10 sequences will be removed;
- With a cutoff of 0.5: a site with 4 gaps over 10 sequences will not be removed;
- With a cutoff of 0.0 a site with 1 gap over 10 sequences will be removed.

#### Usage
```
Usage:
  goalign clean [flags]

Flags:
  -c, --cutoff float    Cutoff for gap deletion : 0 remove sites with > 0 gap, 1 remove sites with 100% gaps)
  -o, --output string   Cleaned alignment output file (default "stdout")
  -q, --quiet           Do not print logs on stderr

Global Flags:
  -i, --align string   Alignment input file (default "stdin")
  -p, --phylip         Alignment is in phylip? False=Fasta
  --input-strict       Strict phylip input format (only used with -p)
  --output-strict      Strict phylip output format  (only used with -p)
```

#### Examples

* Generating a gapped alignment, and removing positions having at least one gap:
```
goalign random -s 10 | goalign mutate gaps -n 1 -r 0.1 -s 10 |  goalign clean
```

Should give on stdout:

```
>Seq0000
ATATGGCGATCAAAGTTCCAATGAGATACTTCCTTTACG
>Seq0001
GCGATTATACTCCGAGGTTTTCCGAATAGTTGTTGGAGT
>Seq0002
TACGTTCGTATAAGCTGCCATCCCGGTGCAGACAATACC
>Seq0003
ATGAGGAAAAATGACGGCTAGCAGCAGGTCGCTCTGTAC
>Seq0004
AACCCCGCCGCACATACGTACTCGCTGTAAGTTCCATTG
>Seq0005
GTGAGCGCTTATGCTTAAGTGATTTGCGTTATATTTCAA
>Seq0006
CCACATTCGTATCTTCCGCCCCCCGGTCCACGCTGCTTC
>Seq0007
TGAACTTTCAGACCAGCGTGCTTTGACACCTGGCTCGGC
>Seq0008
CTAAAACCCTAATCCCTTGCTTCTATCCGCCTGGAGCCT
>Seq0009
CTCGCACGACCAGAGAACTAGGCCGGTTCCGGATATCTA
```

And on stderr:
```
[Warning] message: Alignment (0) length before cleaning=100
[Warning] message: Alignment (0) length after cleaning=39
[Warning] message: Alignment (0) number of gaps=61
```
