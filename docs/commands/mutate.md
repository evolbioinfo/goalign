# Goalign: toolkit and api for alignment manipulation

## Commands

### mutate
This command adds different type of noises in an input alignment, with these sub-commands:
* `goalign mutate gaps` : Adds a given proportion of gaps to a given proportion of the sequences randomly in the input alignment (uniformly).
* `goalign mutate snvs`: Substitute nucleotides/aminoacids by random (uniform) nucleotides/aminoacids with a given rate.

#### Usage
* General command:
```
Usage:
  goalign mutate [command]

Available Commands:
  gaps        Adds gaps uniformly in an input alignment
  snvs        Adds substitutions uniformly in an input alignment

Flags:
  -h, --help            help for mutate
  -o, --output string   Mutated alignment output file (default "stdout")
  -r, --rate float      Mutation rate per nucleotide/amino acid (default 0.1)
  -s, --seed int        Initial Random Seed (default 1497950022161439493)

Global Flags:
  -i, --align string   Alignment input file (default "stdin")
  -p, --phylip         Alignment is in phylip? False=Fasta
```

* gaps command:
```
Usage:
  goalign mutate gaps [flags]

Flags:
  -n, --prop-seq float   Proportion of the sequences in which to add gaps (default 0.5)

Global Flags:
  -i, --align string    Alignment input file (default "stdin")
  -o, --output string   Mutated alignment output file (default "stdout")
  -p, --phylip          Alignment is in phylip? False=Fasta
  -r, --rate float      Mutation rate per nucleotide/amino acid (default 0.1)
  -s, --seed int        Initial Random Seed (default: number of nanoseconds elapsed since January 1, 1970 UTC)
```

* snvs command:
```
Usage:
  goalign mutate snvs [flags]

Global Flags:
  -i, --align string    Alignment input file (default "stdin")
  -o, --output string   Mutated alignment output file (default "stdout")
  -p, --phylip          Alignment is in phylip? False=Fasta
  -r, --rate float      Mutation rate per nucleotide/amino acid (default 0.1)
  -s, --seed int        Initial Random Seed (default: number of nanoseconds elapsed since January 1, 1970 UTC)
```

#### Examples
* Generating a random (uniform) alignment and adding 20% gaps to 50% of the sequences:
```
goalign random -l 20 -s10| goalign mutate gaps -n 0.5 -r 0.2 -s10
```

Should give:
```
>Seq0000
GATTAATTTGCCGTAGGCCA
>Seq0001
G-ATCTGAAGA-CG-A-ACT
>Seq0002
TTAAGTTTT-AC--CTAA-G
>Seq0003
GAGAGGACTAGTTCATACTT
>Seq0004
TT-AAACA-TTTTA-A-CGA
>Seq0005
TGTCGGACCTAAGTATTGAG
>Seq0006
TAC-A-G-TGTATT-CAGCG
>Seq0007
GTGGAGAGGTCTATTTTTCC
>Seq0008
GGTTGAAG-ACT-TA-AGC-
>Seq0009
GTAAAGGGTATGGCCATGTG
```

* Generating a random (uniform) alignment and adding 10% substitutions :
```
goalign random -l 20 -s10| goalign mutate snvs -r 0.1 -s 10
```

Should give:
```
>Seq0000
GATTAATTTCCCGTAGGCCA
>Seq0001
GAATCTGAATATCGAACTAT
>Seq0002
TTAAGTTTTCACTTCTAATG
>Seq0003
GAGAGGACTAGTTCATAATT
>Seq0004
TTTTAACACTTTTACATCGA
>Seq0005
TGTCGGACCTAAGTTTTGTG
>Seq0006
TGCAACGATGTACTCCAGCG
>Seq0007
GTGGAGAGGTCTATTTTTGC
>Seq0008
GGTTAAAGGACTCTATAGCT
>Seq0009
GAAAAGGGTATGGCCATGTG
```
