# Goalign: toolkit and api for alignment manipulation

## Commands

### phase
Although similar to goalign phase, goalign phasent does not take into align translation of input sequences.

This command "phases" input sequences on the basis on either a set of input sequences, or the longest detected orf.
To do so, it will:

1. Search for the longest ORF in the dataset if no reference orf(s) is(are) given;
2. For each sequence of the dataset: will take the sequence in forward and revcomp (if `--reverse` is given),
   align it with all ref orfs, and take the phase (fwd or revcomp) and the reference orf giving the best alignment;
   If no phase gives a good alignment in any reference orf (cutoffs given by `--len-cutoff` and `--match-cutoff`),
   then the sequence flagged as removed;
3. For each sequence, take the Start corresponding to the Start of the ORF, and remove
   nucleotides before (and nucleotides after if `--cut-end` is given);
4. Return the trimmed nucleotidic sequences (phased), the corresponding amino-acid sequences (phasedaa)
   and the start position on the original nt sequence;
5. The log file contains information on:
    1. Sequence name
    2. Its best matching reference orf
    3. Start position on original nt sequence
    4. Extracted sequence length
    5. Positions of nt not in phase with reference orf
    6. Position of the first stop in phase

if `--unaligned` is set, format options are ignored (phylip, nexus, etc.), and
only Fasta is accepted. Otherwise, alignment is first "unaligned".

If input sequences are not nucleotidic, then returns an error.

Output file is an unaligned set of sequences in fasta.

#### Usage
```
Usage:
  goalign phasent [flags]

Flags:
      --aa-output string     Output translated sequences FASTA file (default "none")
      --cut-end              Iftrue, then also remove the end of sequences that do not align with orf
  -h, --help                 help for phasent
      --len-cutoff float     Length cutoff, over orf length, to consider sequence hits (-1==No cutoff) (default -1)
  -l, --log string           Output log: positions of the considered ATG for each sequence (default "none")
      --match float          Score for a match for pairwise alignment (if omitted, then take substitution matrix) (default 1)
      --match-cutoff float   Nb Matches cutoff, over alignment length, to consider sequence hits (-1==No cutoff) (default 0.5)
      --mismatch float       Score for a mismatch for pairwise alignment (if omitted, then take substitution matrix) (default -1)
      --nt-output string     Output ATG "phased" FASTA file + first nts not in ref phase removed (nt corresponding to aa-output sequence) (default "none")
  -o, --output string        Output ATG "phased" FASTA file (default "stdout")
      --ref-orf string       Reference ORF to phase against (if none is given, then will try to get the longest orf in the input data) (default "none")
      --reverse              Search ALSO in the reverse strand (in addition to the forward strand)
      --unaligned            Considers sequences as unaligned and only format fasta is accepted (phylip, nexus,... options are ignored)

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
      --seed int        Random Seed: -1 = nano seconds since 1970/01/01 00:00:00 (default -1)
  -t, --threads int     Number of threads (default 1)
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
goalign phasent -i input.fa --unaligned -o stdout --aa-output aa.fa

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
