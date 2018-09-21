# Goalign: toolkit and api for alignment manipulation

## Commands

### sample
This command samples sites or sequences from an input alignment (fasta by default or phylip with `-p`):
1. `goalign sample sites`: take a random subsequence starting at a random position and with the given length from the input alignment;
2. `goalign sample seqs`: take a random subset of the sequences from an input alignment;
3. `goalign sample rarefy`: Take a new sample taking into accounts weights. Each sequence in the alignment has associated counts. The sum s of the counts represents the number of sequences in the underlying initial dataset. The goal is to downsample (rarefy) the initial dataset, by sampling n sequences from s (n<s), and taking the alignment corresponding to this new sample, i.e by taking only unique (different) sequences from it.

If the input alignment contains several alignments (phylip), will process all of them.

#### Usage
* general command:
```
Usage:
  goalign sample [command]
  
  Available Commands:
    seqs        Samples a subset of sequences from the input alignment
    sites       Take a random subalignment
	  
  Global Flags:
    -i, --align string   Alignment input file (default "stdin")
    -p, --phylip         Alignment is in phylip? False=Fasta
        --input-strict   Strict phylip input format (only used with -p)
        --output-strict  Strict phylip output format  (only used with -p)
```
* seqs command
```
Usage:
  goalign sample seqs [flags]
  
  Flags:
    -n, --nb-seq int      Number of sequences to sample from the alignment (default 1)
    -o, --output string   Sampled alignment output file (default "stdout")
        --seed int        Random Seed: -1 = nano seconds since 1970/01/01 00:00:00 (default -1)
		
  Global Flags:
    -i, --align string   Alignment input file (default "stdin")
    -p, --phylip         Alignment is in phylip? False=Fasta
        --input-strict   Strict phylip input format (only used with -p)
        --output-strict  Strict phylip output format  (only used with -p)

```

* sites command
```
Usage:
  goalign sample sites [flags]
  
  Flags:
    -l, --length int      Length of the random sub alignment (default 10)
    -n, --nsamples int    Number of samples to generate (default 1)
    -o, --output string   Alignment output file (default "stdout")
        --seed int        Random Seed: -1 = nano seconds since 1970/01/01 00:00:00 (default -1)
		  
  Global Flags:
    -i, --align string   Alignment input file (default "stdin")
    -p, --phylip         Alignment is in phylip? False=Fasta
        --input-strict   Strict phylip input format (only used with -p)
        --output-strict  Strict phylip output format  (only used with -p)
```

* rarefy command
```
Usage:
  goalign sample rarefy [flags]

Flags:
  -c, --counts string   Count file (tab separated), one line per sequence: seqname\tcount (default "stdin")
  -n, --nb-seq int      Number of sequences to sample from the repeated dataset (from counts) (default 1)
  -o, --output string   Rarefied alignment output file (default "stdout")
      --seed int        Random Seed: -1 = nano seconds since 1970/01/01 00:00:00 (default -1)

Global Flags:
  -i, --align string    Alignment input file (default "stdin")
      --input-strict    Strict phylip input format (only used with -p)
      --output-strict   Strict phylip output format (only used with -p)
  -p, --phylip          Alignment is in phylip? False=Fasta
```

#### Examples

* Generating a random alignment and taking a subset of the sequences
```
goalign random -l 10 --seed 10 | goalign sample seqs -n 3 --seed 10
```
Should give the following alignment:
```
>Seq0001
CCGTAGGCCA
>Seq0002
GAATCTGAAG
>Seq0008
TTTAAACACT
```

* Generating a random alignment and taking a subsequence from it
```
goalign random -l 10 --seed 10 | goalign sample sites -l 5 --seed 10
```
Should give the following alignment:
```
>Seq0000
TTAAT
>Seq0001
GTAGG
>Seq0002
ATCTG
>Seq0003
CGAAC
>Seq0004
AAGTT
>Seq0005
TTCTA
>Seq0006
GAGGA
>Seq0007
TCATA
>Seq0008
TAAAC
>Seq0009
TACAT
```
