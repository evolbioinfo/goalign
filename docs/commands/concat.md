# Goalign: toolkit and api for alignment manipulation

## Commands

### concat
This command concatenates several alignments in one global alignment. Input alignments may be in phylip or fasta format. If input format is phylip, the file may contain several alignments to concatenate : `goalign concat -i several.phy`. If format is Fasta, all fasta files must be given independently with `goalign concat -i first.fa [second.fa, third.fa, ...]` or `goalign -i none [first.fa, second.Fa, third.fa, ...]`. The order of sequences in alignments may be different, `concat` command will match sequences based on their name.

Examples:
- With a cutoff of 0.5: a site with 5 gaps over 10 sequences will be removed;
- With a cutoff of 0.5: a site with 4 gaps over 10 sequences will not be removed;
- With a cutoff of 0.0 a site with 1 gap over 10 sequences will be removed.

#### Usage
```
Usage:
  goalign concat [flags] [alignment files]

Flags:
  -o, --output string   Alignment output file (default "stdout")

Global Flags:
  -i, --align string   Alignment input file (default "stdin")
```

#### Examples

* Generating a random tree with 5 tips ([Gotree](https://github.com/fredericlemoine/gotree)), simulating 3 alignments from this tree ([seq-gen](https://github.com/rambaut/Seq-Gen), shuffle sequence order, and concatenating them:
```
gotree generate yuletree -l 5 -s 1 -o true_tree.nw
seq-gen -op -mGTR -l500 -z 2 -n 3 true_tree.nw | goalign shuffle seqs -p > alignment.phy
goalign concat -i alignment.phy -p | goalign stats -p
```

It should give the following statistics:
```
length	1500
nseqs	5
avgalleles	1.3220
char	nb	freq
A	1894	0.252533
C	1898	0.253067
G	1788	0.238400
T	1920	0.256000
```
