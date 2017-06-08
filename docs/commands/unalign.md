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
goalign random -s 10 -p | goalign shuffle gaps -s 10 -p
# Gapped alignment => unaligned
goalign random -s 10 -p | goalign shuffle gaps -s 10 -p | goalign unalign -p
```
Should output:
```
  10   100
Seq0000  GATTAATTTG CCGTAGGCCA GAATCTGAAG ATCGAACACT TTAAGTTTTC ACTTCTAATG
Seq0001  --TC------ AA---TTG-- TA-AA-GG-G -ATT-CA--G ------A-G- CT-T-T-T-C
Seq0002  -T-----CGG GC--A-T--T G---GA-CAA GGTT-A---- -----A--GC --CATGATC-
Seq0003  GAGTGGAGGC TTTATGGCAC AAGGTATTAG AGACTGAGGG GCACCCCGGC ATGGTAAGCA
Seq0004  -AT-GC---- -AT--C--G- -C-GT--C-- GG-AA-G--T AC--T-CAC- A---A-AC-C
Seq0005  AGTTTGACTA TGAGCGCCGG CTTAGTGCTG ACAGTGATGC TCCGTTGTAA GGGTCCTGAT
Seq0006  T--C-A---- G-G--A-GT- CGT--T-G-A --AA-CAGCG -----C-CC- ACA------C
Seq0007  CTGGTAATAC CTGCGCTATT TCGTCAGTTC GTGTACGGGT AACGATAGCG GTTAATGCTT
Seq0008  T-G-T--CCC AC--T--C-A CCT------- -G--ATC-G- T-C-C----T GG----CTT-
Seq0009  ACCTACGGCT CTAGACAGCT GAAGTCCGGT TCCGAGCACT GTACGGAAAC TTGAAAAGGC

   GAGAGGACTA GTTCATACTT TTTAAACACT TTTACATCGA
   GG-TGAA--A C----G---T GT-AA---TA -G---ATGTG
   C-C--GG-C- -T--G---GA AG-----A-A -ATGA-CCC-
   GGAGCCATCG CGAAGGCTTC AGGTATCTTC CTGTGTTACC
   CG--GCTA-- CG-CTC---- T-CTT-T-T- C----TCT--
   GTTCTTGTGC TCGCGCATAT TAGAGCTGAG TTTCCCAAAG
   T--GT-GCTC CT--C--C-A ---G-GG-CC TG-G-T--C-
   ATTCCGATCA GCTCACACCC ATGAAGGTGG CTCTGGAGCC
   -T----T-GG -CCCCAGG-- TCA--C--GT -A-CTCT---
   TCGACGGAGG CTTGTTCCGC AGAGTGGGAC TATAACATAC
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
