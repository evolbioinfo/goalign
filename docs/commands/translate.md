# Goalign: toolkit and api for alignment manipulation

## Commands

### translate
Translates an input alignment in amino acids.

If the input alignment is not nucleotides, then returns an error.

It is possible to drop a given number of characters from the start 
of the alignment, by specifying the '--phase' option.

If given phase is -1, then it will translate in the 3 phases, 
from positions 0, 1 and 2. Sequence names will be added the suffix
_<phase>. At the end, 3x times more sequences will be present in the
file.

It only translates using the standard genetic code so far.

#### Usage
```
Usage:
  goalign translate [flags]

Flags:
  -h, --help            help for translate
  -o, --output string   Output translated alignment file (default "stdout")
      --phase int       Number of characters to drop from the start of the alignment (if -1: Translate in the 3 phases, from positions 0, 1, and 2)
      --unaligned       Considers sequences as unaligned and format fasta (phylip, nexus,... options are ignored)

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

* Translating in 3 forward strand phases:

```
goalign translate -i seq.fa --phase -1
```

Should output:
```
>Seq0000_0
FLP
>Seq0000_1
FYP
>Seq0000_2
STP
>Seq0001_0
LKD
>Seq0001_1
LKI
>Seq0001_2
*R*
>Seq0002_0
*HL
>Seq0002_1
NT*
>Seq0002_2
TLE
```
