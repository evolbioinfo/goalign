# GoAlign
GoAlign is a set of command line tools to manipulate phylogenetic trees. It is implemented in [Go](https://golang.org/) language.

## Installation
### Binaries
You can download already compiled binaries for the latest release in the [release](https://github.com/fredericlemoine/goalign/releases) section.
Binaries are available for MacOS, Linux, and Windows (32 and 64 bits).

Once downloaded, you can just run the executable without any other downloads.

### From sources
In order to compile GoAlign, you must first [download](https://golang.org/dl/) and [install](https://golang.org/doc/install) Go on your system.

Then you just have to type :
```
go get github.com/fredericlemoine/goalign/
```
This will download GoAlign sources from github, and all its dependencies.

You can then build it with:
```
cd $GOPATH/src/github.com/fredericlemoine/goalign/
make
```
The `goalign` executable should be located in the `$GOPATH/bin` folder.

## Usage
### List of commands
* addid:      Adds a string to each sequence identifier of the input alignment
* build:       Command to build output files : bootstrap for example
  * seqboot : Generate bootstrap alignments
* clean:       Removes gap sites
* compute:     Different computations (distances, etc.)
  * distances: compute evolutionary distances for nucleotide alignment
* divide:      Divide an input alignment in several output files (one per alignment)
* random:      Generate random sequences
* reformat:    Reformats input alignment into phylip of fasta format
  * fasta
  * nexus
  * phylip
  * tnt
* rename:      Rename sequences of the input alignment, given a map file
* sample:      Randomly samples a subset of sequences from the input alignment
* shuffle:     A set of commands to shuffle an alignment
  * recomb: Recombine some sequences (copy/paste)
  * seqs: Shuffle sequence order in the alignment
  * sites: Shuffle "vertically" some sites of the alignments
  * swap:  Swap portions of some sequences (cut/paste)
* stats:       Prints different characteristics of the alignment
  * alleles
  * char
  * length
  * nalign
  * nseq
  * taxa
* subset:      Take a subset of sequences from the input alignment
* trim:        This command trims names of sequences or sequences themselves
  * name
  * seq
* unalign:     Unaligns input alignment
* version:     Prints the current version of goalign

### Examples

* Generate a random alignemnt and print statistics
```
goalign random | goalign stats
```
* Trim names of a random alignment and finally rename it back
```
goalign random > align.fa
goalign trim name -n 3 -m map -i align.fa > align_rename.fa
goalign rename -i align_rename.fa -m map -r 
```
* Reformat a fasta alignment to phylip
```
goalign random | goalign reformat phylip
```
* Reformat a phylip alignment to fasta
```
goalign random -p | goalign reformat fasta -p
```
* Add a prefix to all sequence names of the alignment
```
goalign random  | goalign addid -n "Dataset1_" 
```
* Add a suffix to all sequence names of the alignment
```
goalign random  | goalign addid -r -n "_Dataset1" 
```
* Take a random sample (10 sequences) from an input alignment
```
goalign random -n 10000 | goalign sample -n 10
```
* Build 100 bootstrap alignments from an input alignment, in a single tar.gz file (5 threads)
```
goalign random -n 500 | goalign build seqboot -S -n 100 --gz --tar -t 5 -o boot
```
* Build 100 bootstrap alignments from an input alignment, in 100 .gz files (5 threads)
```
goalign random -n 500 | goalign build seqboot -S -n 100 --gz -t 5 -o boot
```
