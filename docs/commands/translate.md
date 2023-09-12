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

It is possible to specify alternative genetic code with --genetic-code 
(mitoi, mitov, or standard).

IUPAC codes are taken into account for the translation. If a codon containing 
IUPAC code is ambiguous for translation, then a X is added in place of the aminoacid.

If --ref-seq is given, be careful about the behavior! As with goalign extract, it will will translate
the alignment with the following process: The alignment will be translated codon by
codon using the given reference sequence as guide, by iterating over the reference non gap nucleotides 3 by 3. 
At each iteration, the current reference codon may have gaps between nucleotides, and the translation of the
current codon will be done as following:
	* ex 1:
		Ref: AC--GTACGT
		Seq: ACTTGTACGT
		In that case, the first ref codon is [0,1,4], corresponding to sequence ACTTG in seq
		ACTTG % 3 != 0 ==> Frameshift? => Replaced by T in ref and X in the compared sequence.
	* ex 2:
		Ref: AC---GTACGT
		Seq: ACTTTGTACGT
		ref codon: [0,1,5]
		seq      : ACTTTG (%3==0): Insertion - OK => Replaced by "T-" in ref and "TL" in seq
	* ex 3:
		Ref: ACGTACGT
		Seq: A--TACGT
		ref codon: [0,1,2]
		seq      : A--: Deletion: not ok : Frameshift? => Replaced by "T" in ref and "X" in comp
	* ex 4:
		Ref: AC----GTACGT
		Seq: ACTT-TGTACGT
		ref codon: [0,1,6]
		seq      : ACTTTG (%3==0): Insertion - OK => Replaced by "T-" in ref and "TT" in seq
	* ex 5:
		Ref: AC----GTACGT
		Seq: ACT--TGTACGT
		ref codon: [0,1,6]
		seq      : ACTTTG : Insertion not OK : Frameshift? => Replaced by "T-" in ref and "XX" in seq
This allows to easily translate a multiple sequence alignment containing partial sequences, but the 
interpretation should be careful: the translation of some sequences may not be representative of the 
translation of the unaligned sequences.

#### Usage
```
Usage:
  goalign translate [flags]

Flags:
      --genetic-code string   Genetic Code: standard, mitoi (invertebrate mitochondrial) or mitov (vertebrate mitochondrial) (default "standard")
  -o, --output string         Output translated alignment file (default "stdout")
      --phase int             Number of characters to drop from the start of the alignment (if -1: Translate in the 3 phases, from positions 0, 1, and 2)
      --unaligned             Considers sequences as unaligned and format fasta (phylip, nexus,... options are ignored)

Global Flags:
  -i, --align string          Alignment input file (default "stdin")
      --auto-detect           Auto detects input format (overrides -p, -x and -u)
  -u, --clustal               Alignment is in clustal? default fasta
      --ignore-identical  int Ignore duplicated sequences that have the same name and same sequences
      --input-strict          Strict phylip input format (only used with -p)
  -x, --nexus                 Alignment is in nexus? default fasta
      --no-block              Write Phylip sequences without space separated blocks (only used with -p)
      --one-line              Write Phylip sequences on 1 line (only used with -p)
      --output-strict         Strict phylip output format (only used with -p)
  -p, --phylip                Alignment is in phylip? default fasta
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
