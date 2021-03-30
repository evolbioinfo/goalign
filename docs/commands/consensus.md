# Goalign: toolkit and api for alignment manipulation

## Commands

### consensus
This command generates a basic "majority consensus" sequence, _i.e._ a single sequence whose sites correspond to the majority characters at each positions. 

If '-' is the most abundant character, then '-' will be in the consensus, except if `--ignore-gaps` is specified. If `--ignore-gaps`is specified, then the majority is computed on non gaps characters, except if the column is only made of gaps.
If 'N' is the most abundant character, then 'N' will be in the consensus, except if `--ignore-n` is specified. If `--ignore-n`is specified, then the majority is computed on non N/n characters (X/x for proteins), except if the column is only made of N/n (X/x).


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
      --ignore-gaps        Ignore gaps (except if only gaps on the column)
      --ignore-n           Ignore Ns  (except if only N on the column)
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
