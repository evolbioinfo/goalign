# Goalign: toolkit and api for alignment manipulation

## Commands

### transpose
Transposes an input alignment such that the sequences become the sites and the sites become the sequence.

Each sequence of the output alignment is one site of the input alignment, whose name is the site index (starting from 0).

#### Usage
```
Usage:
  goalign transpose [flags]

Flags:
  -h, --help            help for transpose
  -o, --output string   Output transposed alignment (default "stdout")

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
* Transposing an input alignment:

seq.fa
```
>Seq0000
CTTTC
>Seq0001
GCAAA
>Seq0002
GTGCA
>Seq0003
GTCCG
>Seq0004
TGAGT
```


```
goalign transpose -i seq.fa 
```

Should output:
```
>0
CGGGT
>1
TCTTG
>2
TAGCA
>3
TACCG
>4
CAAGT
```
