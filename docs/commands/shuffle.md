# Goalign: toolkit and api for alignment manipulation

## Commands

### shuffle
This command adds different type of noises in an input alignment, with these sub-commands:
* `goalign shuffle gaps` : Adds a given proportion of gaps to a given proportion of the sequences randomly in the input alignment (uniformly).
* `goalign shuffle recomb`:Recombine a given proportion of the length of a given proportion of the sequences with other sequences (copy/paste). 
    Example:
	```
	s1 CCCCCCCCCCCCCC    s1 CCCCCCCCCCCCCC
	s2 AAAAAAAAAAAAAA => s2 AAAATTTTTTTAAA
	s3 GGGGGGGGGGGGGG    s3 GGGGGGGGGGGGGG
	s4 TTTTTTTTTTTTTT    s4 TTTTTTTTTTTTTT
	```
* `goalign shuffle rogue`: Simulate rogue taxa, by shuffling (horizontally) a given proportion of the sites of a given proportion of the sequences. 
	Example:
	```
	Before        After
	S1 12345678    S1 1234567
	S2 12345678 => S2 1234567
	S3 12345678    S3 1634527 <= This one: 2 nucleotides are shuffled
	S4 12345678    S4 1234567
	```
* `goalign shuffle seqs` : Shuffles sequence order in alignment. It does not add any noise in the alignment.
* `goalign shuffle sites`: Shuffles vertically a given proportion of the alignment sites. Optionnaly, it is possible to shuffle vertically an additional (same) proportion of the remaining intact sites for a given proportion of the sequences (~rogues)
* `goalign shuffle swap` : Swaps a given portion of the sequences in the input alignment. The only option to specify is the swap rate. A rate of 0.5 will swap 25% of the sequences with other 25% of the sequences from a random position to the end. Example swaping 50% of sequences:
  ```
  s1 CCCCCCCCCCCCCC    s1 CCCCCCCCCCCCCC
  s2 AAAAAAAAAAAAAA => s2 AAAAAATTTTTTTT
  s3 GGGGGGGGGGGGGG    s3 GGGGGGGGGGGGGG
  s4 TTTTTTTTTTTTTT    s4 TTTTTTAAAAAAAA
  ```


#### Usage
* General command:
```
Usage:
  goalign shuffle [command]

Available Commands:
  gaps        Adds random (uniformly) gaps
  recomb      Recombine sequences in the input alignment
  rogue       Simulate rogue taxa
  seqs        Shuffles sequence order in alignment
  sites       Shuffles n alignment sites vertically
  swap        Swap portion of sequences in the input alignment

Flags:
  -h, --help            help for shuffle
  -o, --output string   Shuffled alignment output file (default "stdout")
  -s, --seed int        Initial Random Seed (default 1495022409892045891)

Global Flags:
  -i, --align string   Alignment input file (default "stdin")
  -p, --phylip         Alignment is in phylip? False=Fasta
```

* gaps command:
```
Usage:
  goalign shuffle gaps [flags]

Flags:
  -r, --gap-rate float   Proportion of gaps to add per sequences (default 0.5)
  -n, --prop-seq float   Proportion of the sequences to add gaps in (default 0.5)

Global Flags:
  -i, --align string    Alignment input file (default "stdin")
  -o, --output string   Shuffled alignment output file (default "stdout")
  -p, --phylip          Alignment is in phylip? False=Fasta
  -s, --seed int        Initial Random Seed (default: number of nanoseconds elapsed since January 1, 1970 UTC)
```

* recomb command:
```
Usage:
  goalign shuffle recomb [flags]

Flags:
  -l, --prop-length float   Proportion of length of sequences to recombine (default 0.5)
  -n, --prop-seq float      Proportion of the  sequences to recombine (default 0.5)

Global Flags:
  -i, --align string    Alignment input file (default "stdin")
  -o, --output string   Shuffled alignment output file (default "stdout")
  -p, --phylip          Alignment is in phylip? False=Fasta
  -s, --seed int        Initial Random Seed (default 1495021613485404707)
```

* rogue command:
```
Usage:
  goalign shuffle rogue [flags]

Flags:
  -l, --length float        Proportion of the sites to shuffle (default 0.5)
  -n, --prop-seq float      Proportion of the  sequences to consider as rogue (default 0.5)
      --rogue-file string   Rogue sequence names output file (default "stdout")

Global Flags:
  -i, --align string    Alignment input file (default "stdin")
  -o, --output string   Shuffled alignment output file (default "stdout")
  -p, --phylip          Alignment is in phylip? False=Fasta
  -s, --seed int        Initial Random Seed (default 1495022447513391592)
```

* seqs command:
```
Usage:
  goalign shuffle seqs [flags]

Global Flags:
  -i, --align string    Alignment input file (default "stdin")
  -o, --output string   Shuffled alignment output file (default "stdout")
  -p, --phylip          Alignment is in phylip? False=Fasta
  -s, --seed int        Initial Random Seed (default 1495022604538616837)
```

* sites command:
```
Usage:
  goalign shuffle sites [flags]

Flags:
  -r, --rate float          Rate of shuffled sites (>=0 and <=1) (default 0.5)
      --rogue float         If set, then will take the given proportion of taxa, and will apply shuffle again on --rate of the remaining intact sites
      --rogue-file string   Rogue sequence names output file (default "stdout")

Global Flags:
  -i, --align string    Alignment input file (default "stdin")
  -o, --output string   Shuffled alignment output file (default "stdout")
  -p, --phylip          Alignment is in phylip? False=Fasta
  -s, --seed int        Initial Random Seed (default 1495022638306057739)
```

* swap command:
```
Usage:
  goalign shuffle swap [flags]

Flags:
  -r, --rate float   Rate of Swap sequences (>=0 and <=1) (default 0.5)

Global Flags:
  -i, --align string    Alignment input file (default "stdin")
  -o, --output string   Shuffled alignment output file (default "stdout")
  -p, --phylip          Alignment is in phylip? False=Fasta
  -s, --seed int        Initial Random Seed (default 1495022654265868607)
  -t, --threads int     Number of threads (default 1)
```

#### Examples
* Generating a random (uniform) alignment and adding gaps:
```
goalign random -l 20 -s10| goalign shuffle gaps -n 0.5 -r 0.2 -s10
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

* Recombining 0.25% of the sequences on 0.5% of their length
```
# Writing alignment
cat > align.fa <<EOF
>s1
AAAAAAAAAA
>s2
CCCCCCCCCC
>s3
GGGGGGGGGG
>s4
TTTTTTTTTT
EOF
# Recombining alignment
goalign shuffle recomb -i align.fa -l 0.5 -n 0.25 -s 11
```

Should give:
```
>s1
AAAAAAAAAA
>s2
CCCCCCCCCC
>s3
GGGGGGGGGG
>s4
TTCCCCCTTT
  -----
```
