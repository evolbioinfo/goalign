# Goalign: toolkit and api for alignment manipulation

## Commands

### codonalign

Aligns a given nt fasta file using a corresponding aa alignment.

If the input alignment is not amino acid, then returns an error.
If the given fasta file is not nucleotides then returns an error.

Warning: It does not check that the amino acid sequence is a good 
translation of the nucleotide sequence, but just add gaps to the
nucleotide sequence where needed. 

Once gaps are added, if the nucleotide alignment length does not match 
the protein alignment length * 3, returns an error.



#### Usage
```
Usage:
  goalign codonalign [flags]

Flags:
  -f, --fasta string    Input nucleotide Fasta file to be codon aligned (default "stdin")
  -h, --help            help for codonalign
  -o, --output string   Output codon aligned file (default "stdout")

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

input_aa.fa

```
>Seq0000
D*-AVGQNLK
>Seq0001
IE-FKF-LLM
>Seq0002
ERTSSYFLNT
```

input_nt.fa
```
>Seq0000
GATTAAGCCGTAGGCCAGAATCTGAAG
>Seq0001
ATCGAATTTAAGTTTCTTCTAATG
>Seq0002
GAGAGGACTAGTTCATACTTTTTAAACACT
```

```
goalign codonalign -i input_aa.fa -f input_nt.fa 
```

should give

```
>Seq0000
GATTAA---GCCGTAGGCCAGAATCTGAAG
>Seq0001
ATCGAA---TTTAAGTTT---CTTCTAATG
>Seq0002
GAGAGGACTAGTTCATACTTTTTAAACACT
EOF
```
