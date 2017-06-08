# Goalign: toolkit and api for alignment manipulation

## Commands

### unalign
Unaligns an input alignment, by removing indels.

The output is in Fasta format, whatever the input format is (Fasta or Phylip).

Output sequences are free from all indel charachers "-".

As there may be several alignments in the input alignment (phylip format), output files are prefixed with argument given to "--output-prefix", and suffixed with an index and the extension ".fa".

If --output-prefix is set to "-" or "stdout", all sequences are printed on stdout

Example:

goalign unalign -i align.ph -p -o seq_

If align contains 3 alignments, this will generate 3 files:
* seq_000001.fa
* seq_000002.fa
* seq_000003.fa

#### Usage
```
Usage:
  goalign unalign [flags]

Flags:
  -o, --output-prefix string   Unaligned alignment output file prefix (default "stdout")

Global Flags:
  -i, --align string   Alignment input file (default "stdin")
  -p, --phylip         Alignment is in phylip? False=Fasta
  -t, --threads int    Number of threads (default 1)
```

#### Examples

* Generating random alignment, adding 50% gaps in 50% of the sequences, and unaligning it:
```
# Gapped alignment
goalign random -s 10 | goalign shuffle gaps -s 10
# Gapped alignment => unaligned
goalign random -s 10 | goalign shuffle gaps -s 10 | goalign unalign
```
Should output:
```
>Seq0000
GATTAATTTGCCGTAGGCCAGAATCTGAAGATCGAACACTTTAAGTTTTCACTTCTAATGGAGAGGACTAGTTCATACTT
TTTAAACACTTTTACATCGA
>Seq0001
--TC------AA---TTG--TA-AA-GG-G-ATT-CA--G------A-G-CT-T-T-T-CGG-TGAA--AC----G---T
GT-AA---TA-G---ATGTG
>Seq0002
-T-----CGGGC--A-T--TG---GA-CAAGGTT-A---------A--GC--CATGATC-C-C--GG-C--T--G---GA
AG-----A-A-ATGA-CCC-
>Seq0003
GAGTGGAGGCTTTATGGCACAAGGTATTAGAGACTGAGGGGCACCCCGGCATGGTAAGCAGGAGCCATCGCGAAGGCTTC
AGGTATCTTCCTGTGTTACC
>Seq0004
-AT-GC-----AT--C--G--C-GT--C--GG-AA-G--TAC--T-CAC-A---A-AC-CCG--GCTA--CG-CTC----
T-CTT-T-T-C----TCT--
>Seq0005
AGTTTGACTATGAGCGCCGGCTTAGTGCTGACAGTGATGCTCCGTTGTAAGGGTCCTGATGTTCTTGTGCTCGCGCATAT
TAGAGCTGAGTTTCCCAAAG
>Seq0006
T--C-A----G-G--A-GT-CGT--T-G-A--AA-CAGCG-----C-CC-ACA------CT--GT-GCTCCT--C--C-A
---G-GG-CCTG-G-T--C-
>Seq0007
CTGGTAATACCTGCGCTATTTCGTCAGTTCGTGTACGGGTAACGATAGCGGTTAATGCTTATTCCGATCAGCTCACACCC
ATGAAGGTGGCTCTGGAGCC
>Seq0008
T-G-T--CCCAC--T--C-ACCT--------G--ATC-G-T-C-C----TGG----CTT--T----T-GG-CCCCAGG--
TCA--C--GT-A-CTCT---
>Seq0009
ACCTACGGCTCTAGACAGCTGAAGTCCGGTTCCGAGCACTGTACGGAAACTTGAAAAGGCTCGACGGAGGCTTGTTCCGC
AGAGTGGGACTATAACATAC
```

Then:
```
>Seq0000
GATTAATTTGCCGTAGGCCAGAATCTGAAGATCGAACACTTTAAGTTTTCACTTCTAATGGAGAGGACTAGTTCATACTT
TTTAAACACTTTTACATCGA
>Seq0001
TCAATTGTAAAGGGATTCAGAGCTTTTCGGTGAAACGTGTAATAGATGTG
>Seq0002
TCGGGCATTGGACAAGGTTAAGCCATGATCCCGGCTGGAAGAAATGACCC
>Seq0003
GAGTGGAGGCTTTATGGCACAAGGTATTAGAGACTGAGGGGCACCCCGGCATGGTAAGCAGGAGCCATCGCGAAGGCTTC
AGGTATCTTCCTGTGTTACC
>Seq0004
ATGCATCGCGTCGGAAGTACTCACAAACCCGGCTACGCTCTCTTTTCTCT
>Seq0005
AGTTTGACTATGAGCGCCGGCTTAGTGCTGACAGTGATGCTCCGTTGTAAGGGTCCTGATGTTCTTGTGCTCGCGCATAT
TAGAGCTGAGTTTCCCAAAG
>Seq0006
TCAGGAGTCGTTGAAACAGCGCCCACACTGTGCTCCTCCAGGGCCTGGTC
>Seq0007
CTGGTAATACCTGCGCTATTTCGTCAGTTCGTGTACGGGTAACGATAGCGGTTAATGCTTATTCCGATCAGCTCACACCC
ATGAAGGTGGCTCTGGAGCC
>Seq0008
TGTCCCACTCACCTGATCGTCCTGGCTTTTGGCCCCAGGTCACGTACTCT
>Seq0009
ACCTACGGCTCTAGACAGCTGAAGTCCGGTTCCGAGCACTGTACGGAAACTTGAAAAGGCTCGACGGAGGCTTGTTCCGC
AGAGTGGGACTATAACATAC
```
