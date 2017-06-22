# Goalign: toolkit and api for alignment manipulation

## Commands

### addid
This command adds an indentifier (string) to all sequences of an input alignment. The string may be added to the left or to the right of each sequence name.

By default the string is added to the left of each name.

#### Usage

General command
```
Usage:
  goalign addid [flags]

Flags:
  -n, --name string        String to add to sequence names (default "none")
  -o, --out-align string   Renamed alignment output file (default "stdout")
  -r, --right              Adds the String on the right of sequence names (otherwise, adds to left)

Global Flags:
  -i, --align string   Alignment input file (default "stdin")
  -p, --phylip         Alignment is in phylip? False=Fasta
  --input-strict       Strict phylip input format (only used with -p)
  --output-strict      Strict phylip output format  (only used with -p)
```

#### Examples

* Generate a random alignment with 5 sequences, adding "prefix_" as prefix and "_suffix" as suffix to each sequence name.

```
goalign random -s 10 -n 5 | goalign addid -n prefix_ | goalign addid -n _suffix -r
```

Should give

```
>prefix_Seq0000_suffix
GATTAATTTGCCGTAGGCCAGAATCTGAAGATCGAACACTTTAAGTTTTCACTTCTAATGGAGAGGACTAGTTCATACTT
TTTAAACACTTTTACATCGA
>prefix_Seq0001_suffix
TGTCGGACCTAAGTATTGAGTACAACGGTGTATTCCAGCGGTGGAGAGGTCTATTTTTCCGGTTGAAGGACTCTAGAGCT
GTAAAGGGTATGGCCATGTG
>prefix_Seq0002_suffix
CTAAGCGCGGGCGGATTGCTGTTGGAGCAAGGTTAAATACTCGGCAATGCCCCATGATCCCCCAAGGACAATAAGAGCGA
AGTTAGAACAAATGAACCCC
>prefix_Seq0003_suffix
GAGTGGAGGCTTTATGGCACAAGGTATTAGAGACTGAGGGGCACCCCGGCATGGTAAGCAGGAGCCATCGCGAAGGCTTC
AGGTATCTTCCTGTGTTACC
>prefix_Seq0004_suffix
CATAGCCCCTGATGCCCTGACCCGTGTCGCGGCAACGTCTACATTTCACGATAAATACTCCGCTGCTAGTCGGCTCTAGA
TGCTTTTCTTCCAGATCTGG
```
