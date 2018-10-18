# Goalign: toolkit and api for alignment manipulation

## Commands

### phase
This command "phases" sequences on the ATG giving the longest ORF on the forward strand. To do so, phase will:

1. Search for the longest ORF in the dataset if no reference orf is given;
1. Translate the given ORF in aminoacids;
2. For each sequence of the dataset: translate it in the 3 phases (forward strand only),
   align it with the translated orf, and take the phase giving the best alignment; If no phase
   gives a good alignment (>80% orf length, and starting at first position of the ORF), then
   the sequence is discarded;
3. For each sequence, take the Start corresponding to the Start of the ORF, and remove
   nucleotides before;
4. Return the trimmed nucleotidic sequences (phased), the corresponding amino-acid sequences (phasedaa)
   and the positions of starts in the nucleotidic sequences.


if --unaligned is set, format options are ignored (phylip, nexus, etc.), and
only Fasta is accepted.

If input sequences are not nucleotidic, then returns an error.

Output file is an unaligned set of sequences in fasta.

#### Usage
```
Usage:
  goalign phase [flags]

Flags:
      --aa-output string   Output Met "phased" aa FASTA file (default "none")
  -h, --help               help for phase
  -l, --log string         Output log: positions of the considered ATG for each sequence (default "none")
  -o, --output string      Output ATG "phased" FASTA file (default "stdout")
      --ref-orf string     Reference ORF to phase against (if none is given, then will try to get the longest orf in the input data) (default "none")
      --unaligned          Considers sequences as unaligned and only fasta format is accepted (phylip, nexus,... options are ignored)

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

input.fa

```
>allcodons
GCTGCCGCAGCGTTATTGCTTCTCCTACTGCGTCGCCGACGGAGAAGGAAAAAGAATAACATG
GATGACTTTTTCTGTTGCCCTCCCCCACCGCAACAGTCTTCCTCATCGAGTAGCGAAGAGACT
ACCACAACGGGTGGCGGAGGGTGGCATCACTATTACATTATCATAGTTGTCGTAGTGTAATGA
TAGC
>allcodons2
GCTGCAGCGTTATTGCTTCTCCTACTGCGTCGCCGACGGAGAAGGAAAAAGAATAACATG
GATGACTTTTTCTGTTGCCCTCCCCCACCGCAACAGTCTTCCTCATCGAGTAGCGAAGAGACT
ACCACAACGGGTGGCGGAGGGTGGCATCACTATTACATTATCATAGTTGTCGTATAATGA
```

```
goalign phase -i input.fa --unaligned -o stdout --aa-output aa.fa

```

should give

stdout
```
>allcodons
ATGGATGACTTTTTCTGTTGCCCTCCCCCACCGCAACAGTCTTCCTCATCGAGTAGCGAAGAGACTACCACAACGGGTGG
CGGAGGGTGGCATCACTATTACATTATCATAGTTGTCGTAGTGTAATGATAGC
>allcodons2
ATGGATGACTTTTTCTGTTGCCCTCCCCCACCGCAACAGTCTTCCTCATCGAGTAGCGAAGAGACTACCACAACGGGTGG
CGGAGGGTGGCATCACTATTACATTATCATAGTTGTCGTATAATGA
```

aa.fa
```
>allcodons
MDDFFCCPPPPQQSSSSSSEETTTTGGGGWHHYYIIIVVVV***
>allcodons2
MDDFFCCPPPPQQSSSSSSEETTTTGGGGWHHYYIIIVVV**
```
