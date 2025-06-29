# Goalign: toolkit and api for alignment manipulation

## Commands

### mutate
This command adds different type of noises in an input alignment, with these sub-commands:
* `goalign mutate gaps` : Adds a given proportion of gaps to a given proportion of the sequences randomly in the input alignment (uniformly).
* `goalign mutate snvs`: Substitute nucleotides/aminoacids by random (uniform) nucleotides/aminoacids with a given rate. Does not apply to gaps or other special characters.
* `goalign mutate ambig`: Adds a given proportion of ambiguities (N or X depending on the alphabet) to a given proportion of the input sequences.

#### Usage
* General command:
```
Usage:
  goalign mutate [command]

Available Commands:
  ambig       Adds ambiguities uniformly to an input alignment
  gaps        Adds gaps uniformly in an input alignment
  snvs        Adds substitutions uniformly in an input alignment

Flags:
  -h, --help            help for mutate
  -o, --output string   Mutated alignment output file (default "stdout")
  -r, --rate float      Mutation rate per nucleotide/amino acid (default 0.1)

Global Flags:
  -i, --align string           Alignment input file (default "stdin")
      --alphabet string        Alignment/Sequences alphabet: auto (default), aa, or nt (default "auto")
      --auto-detect            Auto detects input format (overrides -p, -x and -u)
  -u, --clustal                Alignment is in clustal? default fasta
      --ignore-identical int   Ignore duplicated sequences that have the same name and potentially have same sequences, 0 : Does not ignore anything, 1: Ignore sequences having the same name (keep the first one whatever their sequence), 2: Ignore sequences having the same name and the same sequence
      --input-strict           Strict phylip input format (only used with -p)
  -x, --nexus                  Alignment is in nexus? default fasta
      --no-block               Write Phylip sequences without space separated blocks (only used with -p)
      --one-line               Write Phylip sequences on 1 line (only used with -p)
      --output-strict          Strict phylip output format (only used with -p)
  -p, --phylip                 Alignment is in phylip? default fasta
      --seed int               Random Seed: -1 = nano seconds since 1970/01/01 00:00:00 (default -1)
  -k, --stockholm              Alignment is in stockholm? default fasta
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
  -   --seed int        Random Seed: -1 = nano seconds since 1970/01/01 00:00:00 (default -1)
  --input-strict        Strict phylip input format (only used with -p)
  --output-strict       Strict phylip output format  (only used with -p)
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
      --seed int        Random Seed: -1 = nano seconds since 1970/01/01 00:00:00 (default -1)
  --input-strict        Strict phylip input format (only used with -p)
  --output-strict       Strict phylip output format  (only used with -p)
```

* ambig command:
```
Usage:
  goalign mutate ambig [flags]

Flags:
  -h, --help             help for ambig
  -n, --prop-seq float   Proportion of the sequences in which to add ambiguities (default 0.5)

Global Flags:
  -i, --align string           Alignment input file (default "stdin")
      --alphabet string        Alignment/Sequences alphabet: auto (default), aa, or nt (default "auto")
      --auto-detect            Auto detects input format (overrides -p, -x and -u)
  -u, --clustal                Alignment is in clustal? default fasta
      --ignore-identical int   Ignore duplicated sequences that have the same name and potentially have same sequences, 0 : Does not ignore anything, 1: Ignore sequences having the same name (keep the first one whatever their sequence), 2: Ignore sequences having the same name and the same sequence
      --input-strict           Strict phylip input format (only used with -p)
  -x, --nexus                  Alignment is in nexus? default fasta
      --no-block               Write Phylip sequences without space separated blocks (only used with -p)
      --one-line               Write Phylip sequences on 1 line (only used with -p)
  -o, --output string          Mutated alignment output file (default "stdout")
      --output-strict          Strict phylip output format (only used with -p)
  -p, --phylip                 Alignment is in phylip? default fasta
  -r, --rate float             Mutation rate per nucleotide/amino acid (default 0.1)
      --seed int               Random Seed: -1 = nano seconds since 1970/01/01 00:00:00 (default -1)
  -k, --stockholm              Alignment is in stockholm? default fasta
```


#### Examples
* Generating a random (uniform) alignment and adding 20% gaps to 50% of the sequences:
```
goalign random -l 20 --seed 10| goalign mutate gaps -n 0.5 -r 0.2 --seed 10
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

* Generating a random (uniform) nucleotide alignment and adding 10% ambiguities :
```
goalign random -l 20 --seed 10| goalign mutate ambig -r 0.8 --seed 10
```

Should give:
```
>Seq0000
GATTAATTTGCCGTAGGCCA
>Seq0001
NNNTNNNNANANNNNNNANN
>Seq0002
NNNNGNNNNNANNNNTNANN
>Seq0003
GAGAGGACTAGTTCATACTT
>Seq0004
TNNNNNNANNNNNNNNNCNA
>Seq0005
TGTCGGACCTAAGTATTGAG
>Seq0006
NNNNNNNNNNTNTTNNNNCN
>Seq0007
GTGGAGAGGTCTATTTTTCC
>Seq0008
GGNNNNNNNNNNNTANNNNN
>Seq0009
GTAAAGGGTATGGCCATGTG
```


* Generating a random (uniform) amino-acid alignment and adding 10% ambiguities :

```
$ goalign random -a -l 20 --seed 10| .goalign mutate ambig -r 0.8 --seed 10
>Seq0000
PHGVHCVSSYRFEKCPNFFC
>Seq0001
XXXKXXXXCXMXXXXXXHXX
>Seq0002
XXXXEXXXXXAXXXXGXHXX
>Seq0003
YHPTYLHWSAPDGRCKTQSV
>Seq0004
DXXXXXXMXXXXXXXXXQXH
>Seq0005
GLKQYYAQRKATNSHKDLAY
>Seq0006
XXXXXXXXXXKXVGXXXXWX
>Seq0007
NSNLHNANYSQGHVVSDVIF
>Seq0008
YYXXXXXXXXXXXSCXXXXX
>Seq0009
PGTTAYLLGHDYNWFCSEKN
```
