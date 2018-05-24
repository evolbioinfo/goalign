# Goalign: toolkit and api for alignment manipulation

## Commands

### translate
This command translates an input sequence into amino acids.


#### Usage
```
Usage:
  goalign translate [flags]

Flags:
  -h, --help            help for translate
  -o, --output string   Output translated alignment file (default "stdout")
      --phase int       Number of characters to drop from the start of the alignment

Global Flags:
  -i, --align string    Alignment input file (default "stdin")
      --auto-detect     Auto detects input format (overrides -p and -x)
      --input-strict    Strict phylip input format (only used with -p)
  -x, --nexus           Alignment is in nexus? default fasta
      --output-strict   Strict phylip output format (only used with -p)
  -p, --phylip          Alignment is in phylip? default fasta
  -t, --threads int     Number of threads (default 1)
```


#### Examples
* Translating an input sequence:

seq.fa
```
>Seq0000
TTTCTACCCCA
>Seq0001
CTTAAAGATAG
>Seq0002
TAACACTTGAA
```


```
goalign translate -i seq.fa --phase 1
```

Should output:
```
>Seq0000
FYP
>Seq0001
LKI
>Seq0002
NT*
```
