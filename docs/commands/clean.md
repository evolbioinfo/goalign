# Goalign: toolkit and api for alignment manipulation

## Commands

### clean
This command removes alignment sites or sequences constitued of >= than a given proportion of gaps. Exception for a cutoff of 0: removes sites /sequences constitued of > 0 gaps.

Two subcommands:

* `goalign clean sites`: Removes sites with gaps
* `goalign clean seqs`: Removes sequences with gaps

Examples with sites:
- With a cutoff of 0.5: a site with 5 gaps over 10 sequences will be removed;
- With a cutoff of 0.5: a site with 4 gaps over 10 sequences will not be removed;
- With a cutoff of 0.0 a site with 1 gap over 10 sequences will be removed.

#### Usage
```
Usage:
  goalign clean [command]

Available Commands:
  seqs        Removes sequences with gaps
  sites       Removes sites with gaps

Flags:
  -c, --cutoff float    Cutoff for gap deletion : 0 remove sites/sequences with > 0 gap, 1 remove sites/sequences with 100% gaps)
  -h, --help            help for clean
  -o, --output string   Cleaned alignment output file (default "stdout")
  -q, --quiet           Do not print results on stderr

Global Flags:
  -i, --align string    Alignment input file (default "stdin")
      --input-strict    Strict phylip input format (only used with -p)
      --output-strict   Strict phylip output format (only used with -p)
  -p, --phylip          Alignment is in phylip? False=Fasta
```

#### Examples

* Generating a gapped alignment, and removing positions having at least one gap:
```
goalign random --seed 10 | goalign mutate gaps -n 1 -r 0.1 --seed 10 |  goalign clean sites
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

* Generating a gapped alignment, and removing sequences having at least one gap:
```
goalign random --seed 10 | goalign mutate gaps -n 0.5 -r 0.7 --seed 10 |  goalign clean seqs
```

Should give on stdout:

```
>Seq0000
GATTAATTTGCCGTAGGCCAGAATCTGAAGATCGAACACTTTAAGTTTTCACTTCTAATGGAGAGGACTAGTTCATACTT
TTTAAACACTTTTACATCGA
>Seq0003
GAGTGGAGGCTTTATGGCACAAGGTATTAGAGACTGAGGGGCACCCCGGCATGGTAAGCAGGAGCCATCGCGAAGGCTTC
AGGTATCTTCCTGTGTTACC
>Seq0005
AGTTTGACTATGAGCGCCGGCTTAGTGCTGACAGTGATGCTCCGTTGTAAGGGTCCTGATGTTCTTGTGCTCGCGCATAT
TAGAGCTGAGTTTCCCAAAG
>Seq0007
CTGGTAATACCTGCGCTATTTCGTCAGTTCGTGTACGGGTAACGATAGCGGTTAATGCTTATTCCGATCAGCTCACACCC
ATGAAGGTGGCTCTGGAGCC
>Seq0009
ACCTACGGCTCTAGACAGCTGAAGTCCGGTTCCGAGCACTGTACGGAAACTTGAAAAGGCTCGACGGAGGCTTGTTCCGC
AGAGTGGGACTATAACATAC
```

And on stderr:
```
[Warning] in cmd/cleanseqs.go (line 36), message: Alignment (0) #seqs before cleaning=10
[Warning] in cmd/cleanseqs.go (line 37), message: Alignment (0) #seqs after cleaning=5
[Warning] in cmd/cleanseqs.go (line 38), message: Alignment (0) removed sequences=5
```
