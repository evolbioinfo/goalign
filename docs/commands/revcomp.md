# Goalign: toolkit and api for alignment manipulation

## Commands

### revcomp
Reverse complements an input alignment.

If the input alignment is not nucleotides, then returns an error.

If `--unaligned` is specified, then input sequences may be unaligned.

IUPAC codes are taken into account.

#### Usage
```
Usage:
  goalign revcomp [flags]

Flags:
  -h, --help            help for revcomp
  -o, --output string   Output reverse complement alignment file (default "stdout")
      --unaligned       Considers sequences as unaligned and format fasta (phylip, nexus,... options are ignored)

Global Flags:
  -i, --align string          Alignment input file (default "stdin")
      --auto-detect           Auto detects input format (overrides -p, -x and -u)
  -u, --clustal               Alignment is in clustal? default fasta
      --ignore-identical int  Ignore duplicated sequences that have the same name and same sequences
      --input-strict          Strict phylip input format (only used with -p)
  -x, --nexus                 Alignment is in nexus? default fasta
      --no-block              Write Phylip sequences without space separated blocks (only used with -p)
      --one-line              Write Phylip sequences on 1 line (only used with -p)
      --output-strict         Strict phylip output format (only used with -p)
  -p, --phylip                Alignment is in phylip? default fasta
```


#### Examples
* Reverse complement:

seq.fa
```
>Seq0000
CTTTCGCAAA
>Seq0001
GTGCAGTCCG
>Seq0002
TGAGTTTAGT
>Seq0003
CATTCACTCG
>Seq0004
CGGTCTGATC
>Seq0005
CCCTACAGTT
>Seq0006
TGCAGACGTG
>Seq0007
TAGGTGCTAA
>Seq0008
TCCCCTCTTG
>Seq0009
GAGTATATCG
```


```
goalign revcomp -i seq.fa
```

Should output:
```
>Seq0000
TTTGCGAAAG
>Seq0001
CGGACTGCAC
>Seq0002
ACTAAACTCA
>Seq0003
CGAGTGAATG
>Seq0004
GATCAGACCG
>Seq0005
AACTGTAGGG
>Seq0006
CACGTCTGCA
>Seq0007
TTAGCACCTA
>Seq0008
CAAGAGGGGA
>Seq0009
CGATATACTC
```
