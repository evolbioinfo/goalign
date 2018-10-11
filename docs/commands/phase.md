# Goalign: toolkit and api for alignment manipulation

## Commands

### phase
This command "phases" sequences on the ATG giving the longest ORF on the forward strand.

if --unaligned is set, format options are ignored (phylip, nexus, etc.), and
only Fasta is accepted.

If input sequences are not nucleotidic, then returns an error.

Output file is an unaligned set of sequences in fasta.

#### Usage
```
Usage:
  goalign phase [flags]

Flags:
  -l, --log string      Output log: positions of the considered ATG for each sequence (default "none")
  -o, --output string   Output ATG "phased" file (default "stdout")
      --unaligned       Considers sequences as unaligned and only format fasta is accepted (phylip, nexus,... options are ignored)

Global Flags:
  -i, --align string    Alignment input file (default "stdin")
      --auto-detect     Auto detects input format (overrides -p, -x and -u)
  -u, --clustal         Alignment is in clustal? default fasta
      --input-strict    Strict phylip input format (only used with -p)
  -x, --nexus           Alignment is in nexus? default fasta
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
goalign phase -i input --unaligned
```

should give

```
>allcodons
ATGGATGACTTTTTCTGTTGCCCTCCCCCACCGCAACAGTCTTCCTCATCGAGTAGCGAAGAGACTACCACAACGGGTGG
CGGAGGGTGGCATCACTATTACATTATCATAGTTGTCGTAGTGTAATGATAGC
>allcodons2
ATGGATGACTTTTTCTGTTGCCCTCCCCCACCGCAACAGTCTTCCTCATCGAGTAGCGAAGAGACTACCACAACGGGTGG
CGGAGGGTGGCATCACTATTACATTATCATAGTTGTCGTATAATGA
```
