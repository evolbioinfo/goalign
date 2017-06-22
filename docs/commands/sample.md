# Goalign: toolkit and api for alignment manipulation

## Commands

### sample
This command samples sites or sequences from an input alignment (fasta by default or phylip with `-p`):
1. `goalign sample sites`: take a random subsequence starting at a random position and with the given length from the input alignment;
2. `goalign sample seqs`: take a random subset of the sequences from an input alignment;

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
      --input-strict       Strict phylip input format (only used with -p)
      --output-strict      Strict phylip output format  (only used with -p)
```
* seqs command
```
Usage:
  goalign sample seqs [flags]
  
  Flags:
    -n, --nb-seq int      Number of sequences to sample from the alignment (default 1)
	-o, --output string   Sampled alignment output file (default "stdout")
	-s, --seed int        Initial Random Seed (default 1496330366914852296)
		
  Global Flags:
	-i, --align string   Alignment input file (default "stdin")
    -p, --phylip         Alignment is in phylip? False=Fasta
    --input-strict       Strict phylip input format (only used with -p)
    --output-strict      Strict phylip output format  (only used with -p)

```

* sites command
```
Usage:
  goalign sample sites [flags]
  
  Flags:
    -l, --length int      Length of the random sub alignment (default 10)
	-n, --nsamples int    Number of samples to generate (default 1)
	-o, --output string   Alignment output file (default "stdout")
	-s, --seed int        Initial Random Seed (default 1496330424289100016)
		  
  Global Flags:
    -i, --align string   Alignment input file (default "stdin")
	-p, --phylip         Alignment is in phylip? False=Fasta
    --input-strict       Strict phylip input format (only used with -p)
    --output-strict      Strict phylip output format  (only used with -p)
```

#### Examples

* Generating a random alignment and taking a subset of the sequences
```
goalign random -l 10 -s 10 | goalign sample seqs -n 3 -s 10
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
goalign random -l 10 -s 10 | goalign sample sites -l 5 -s 10
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
