# Goalign: toolkit and api for alignment manipulation

## Commands

### divide
This command divides an input alignment file (phylip format) containing several alignments, in multiple files containing one alignment. If alignment is in Fasta, will create only one file. Option `-o` indicates prefix of output files. 

#### Usage
```
Usage:
  goalign divide [flags]

Flags:
  -f, --out-fasta       Output files in fasta format (default, same as input)
  -o, --output string   Divided alignment output files prefix (default "prefix")

Global Flags:
  -i, --align string   Alignment input file (default "stdin")
  -p, --phylip         Alignment is in phylip? False=Fasta
```

#### Examples

* Generating a random tree with 5 tips ([Gotree](https://github.com/fredericlemoine/gotree)), simulating 3 alignments from this tree ([seq-gen](https://github.com/rambaut/Seq-Gen)), and writing them in 3 independ fasta files:
```
gotree generate yuletree -l 5 -s 1 -o true_tree.nw
seq-gen -op -mGTR -l500 -z 2 -n 3 true_tree.nw | goalign divide -p -f -o align
```

There should be three files in the current directory:
* `align_000.fa`
* `align_001.fa`
* `align_002.fa`
