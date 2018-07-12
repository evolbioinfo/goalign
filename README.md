# Goalign
[![build](https://travis-ci.org/fredericlemoine/goalign.svg?branch=master)](https://travis-ci.org/fredericlemoine/goalign) [![Anaconda-Server Badge](https://anaconda.org/bioconda/goalign/badges/installer/conda.svg)](https://anaconda.org/bioconda/goalign) [![Docker hub](https://img.shields.io/docker/build/evolbioinfo/goalign.svg)](https://hub.docker.com/r/evolbioinfo/goalign/builds/)


Goalign is a set of command line tools to manipulate multiple alignments. It is implemented in [Go](https://golang.org/) language.

Goalign aims to handle multiple alignments in [Phylip](https://en.wikipedia.org/wiki/PHYLIP), [Fasta](https://en.wikipedia.org/wiki/FASTA_format), [Nexus](https://en.wikipedia.org/wiki/Nexus_file), and [Clustal](https://en.wikipedia.org/wiki/Clustal) formats, through several basic commands. Each command may print result (an alignment for example) in the standard output, and thus can be piped to the standard input of the next goalign command.

Input files may be local or remote files:

- If file name is of the form `http(s)://<URL>`, the file is download from the given URL.
- Otherwise, the file is considered local.

Gzipped input files (`.gz` extension) are supported.


**Note**:

TO manipulate phylogenetic trees, See also [Gotree](https://github.com/fredericlemoine/gotree).

## Installation
### Easy way: Binaries
You can download ready to run binaries for the latest release in the [release](https://github.com/fredericlemoine/goalign/releases) section.
Binaries are available for MacOS, Linux, and Windows (32 and 64 bits).

Once downloaded, you can just run the executable without any other downloads.

### Docker
Goalign Docker image is accessible from [docker hub](https://hub.docker.com/r/evolbioinfo/goalign/). You may use it as following:

```[bash]
# Display goalign help
docker run -v $PWD:$PWD -w $PWD -i -t evolbioinfo/goalign:v0.2.6 -h
```

### Singularity
Goalign [docker image](https://hub.docker.com/r/evolbioinfo/goalign/) is usable from singularity . You may use it as following:

```[bash]
# Pull image from docker hub
singularity pull docker://evolbioinfo/goalign:v0.2.6
# Display goalign help
./goalign-v0.2.6.simg -h
```

### Conda
Goalign is also available on [bioconda](https://anaconda.org/bioconda/goalign). Just type:

```
conda install -c bioconda goalign
```

### From sources
In order to compile Goalign, you must first [download](https://golang.org/dl/) and [install](https://golang.org/doc/install) Go on your system.

Then you just have to type :
```
go get github.com/fredericlemoine/goalign/
go get -u github.com/golang/dep/cmd/dep
```
This will download GoAlign sources from github, and all its dependencies.

You can then build it with:
```
cd $GOPATH/src/github.com/fredericlemoine/goalign/
dep ensure
make
```
The `goalign` executable should be located in the `$GOPATH/bin` folder.

## Auto completion

### Bash
* Install bash-completion:
```
# MacOS brew
brew install bash-completion
# MacOS port
port install bash-completion
# Linux /centos
yum install bash-completion -y
# Lunix /ubuntu
apt-get install bash-completion
```

* Activate goalign bash completion
```
# Once
source <(goalign completion bash)
# Permanently
mkdir ~/.goalign
goalign completion bash > ~/.goalign/completion.bash.inc
printf "
# goalign shell completion
source '$HOME/.goalign/completion.bash.inc'
" >> $HOME/.bashrc
```

### Zsh (not tested)

```
# Once
source <(kubectl completion zsh)
# Permanently
goalign completion zsh > "${fpath[1]}/_goalign"
```


## Usage
You may go to the [doc](docs/index.md) for a more detailed documentation of the commands.

### List of commands
* addid:      Adds a string to each sequence identifier of the input alignment
* build:       Command to build output files : bootstrap for example
  * seqboot : Generate bootstrap alignments
* clean:       Removes gap sites/sequences
  * sites : Removes sites with gaps
  * seqs : Removes sequences with gaps
* compute:     Different computations (distances, etc.)
  * distances: compute evolutionary distances for nucleotide alignment
  * entropy: compute entropy of alignment sites
  * pssm: compute position-specific scoring matrix
* concat:      Concatenates several alignments
* divide:      Divide an input alignment in several output files (one per alignment)
* draw:   Draw alignments
  * biojs:     Display an input alignment in an html file using [BioJS](http://msa.biojs.net/)
* mutate: Adds substitutions (~sequencing errors), or gaps, uniformly in an input alignment
  * gaps: Adds gaps uniformly in an input alignment
  * snvs: Adds substitutions uniformly in an input alignment
* random:      Generate random sequences
* reformat:    Reformats input alignment into several formats
  * fasta
  * nexus
  * paml
  * clustal
  * phylip
  * tnt
* rename:      Rename sequences of the input alignment, given a map file
* sample: Samples sequences or subalignments
  * seqs: Randomly samples a subset of sequences from the input alignment
  * sites: Extracts a sub-alignment starting a a random position, and with a given length
  * rarefy: Down-samples input alignment, taking into accounts weights/counts of all sequences
* shuffle:     A set of commands to shuffle an alignment
  * recomb: Recombine some sequences (copy/paste)
  * rogue: simulate sort of rogue taxa by shuffling some sequences
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
* subseq:      Extract a subsequence from the alignment
* subset:      Take a subset of sequences from the input alignment
* trim:        This command trims names of sequences or sequences themselves
  * name
  * seq
* unalign:     Unaligns input alignment
* version:     Prints the current version of goalign

### Goalign commandline examples

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
* Reformat a clustal alignment to fasta
```
goalign random --amino-acids --clustal --nb-seqs 2 | goalign reformat fasta --clustal
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
* Extract all sequences whose name starts with "mammal"
```
goalign subset -e '^mammal.*$' -i align.fasta
```
* Extract all sequences whose name does match the regexp
```
goalign subset -r -e '^mammal.*$' -i align.fasta
```
* Extract a sub sequences going from position 10 and with a length of 100
```
goalign subseq -i align.fasta -s 9 -l 10
```
* Compute a "logo" like consensus
```
goalign compute pssm -n 4 -i align.fasta
```
* Compute an evolutionary distance matrix (dna alignment only, 5 threads)
```
goalign compute distance -m k2p -i align.fasta -t 5
```
* Compute site entropry
```
goalign compute entropy -i align.fasta
```
* Build 100 bootstrap alignments from an input alignment, in a single tar.gz file (5 threads)
```
goalign random -n 500 | goalign build seqboot -S -n 100 --gz --tar -t 5 -o boot
```
* Build 100 bootstrap alignments from an input alignment, in 100 .gz files (5 threads)
```
goalign random -n 500 | goalign build seqboot -S -n 100 --gz -t 5 -o boot
```


### Goalign api usage examples
* Parse a Phylip single alignment file and export it in Fasta
```go
package main

import (
	"fmt"
	"os"

	"github.com/fredericlemoine/goalign/align"
	"github.com/fredericlemoine/goalign/io/fasta"
	"github.com/fredericlemoine/goalign/io/phylip"
)

func main() {
	var err error
	var f *os.File
	var align align.Alignment

	f, err = os.Open("f.phy")
	if err != nil {
		panic(err)
	}
	if align, err = phylip.NewParser(f).Parse(); err != nil {
		panic(err)
	} else {
		fmt.Println(fasta.WriteSequences(align))
	}
}
```

* Parse a Phylip multi alignments file and export it in Fasta
```go
package main

import (
	"fmt"
	"os"

	"github.com/fredericlemoine/goalign/align"
	"github.com/fredericlemoine/goalign/io/fasta"
	"github.com/fredericlemoine/goalign/io/phylip"
)

func main() {
	var f *os.File
	var aligns chan align.Alignment
	var err error

	f, err = os.Open("f.phy")
	if err != nil {
		panic(err)
	}
	aligns = make(chan align.Alignment, 15)
	if err = phylip.NewParser(f).ParseMultiple(aligns); err != nil {
		panic(err)
	} else {
		for al := range aligns {
			fmt.Println(fasta.WriteSequences(al))
		}
	}
}
```
* Parse a Fasta file and export it in Nexus
```go
package main

import (
	"fmt"
	"os"

	"github.com/fredericlemoine/goalign/align"
	"github.com/fredericlemoine/goalign/io/fasta"
	"github.com/fredericlemoine/goalign/io/nexus"
)

func main() {
	var f *os.File
	var align align.Alignment
	var err error

	f, err = os.Open("f.fasta")
	if err != nil {
		panic(err)
	}
	if align, err = fasta.NewParser(f).Parse(); err != nil {
		panic(err)
	} else {
		fmt.Println(nexus.WriteAlignment(align))
	}
}
```
* Parse a Fasta file and export it in Phylip
```go
package main

import (
	"fmt"
	"os"

	"github.com/fredericlemoine/goalign/align"
	"github.com/fredericlemoine/goalign/io/fasta"
	"github.com/fredericlemoine/goalign/io/phylip"
)

func main() {
	var f *os.File
	var align align.Alignment
	var err error

	f, err = os.Open("f.fasta")
	if err != nil {
		panic(err)
	}
	if align, err = fasta.NewParser(f).Parse(); err != nil {
		panic(err)
	} else {
		fmt.Println(phylip.WriteAlignment(align, false))
	}
}
```
* Iterating over alignment sequences
```go
	align.IterateChar(func(name string, sequence []rune) {
		fmt.Printf("Sequence: %s\n", name)
	})
```
* Append identifier at the beginning of all sequence names
```go
align.AppendSeqIdentifier("IDENT", false)
```
* Alignment statistics
```go
var n int = align.NbSequences()
var l int = align.Length()
```
* Extract a sub alignment
```go
var subalign align.Alignment
var err error
subalign,err = align.SubAlign(0, 100)
```
* Sort sequences by alphanumerical order
```go
align.Sort()
```
* Copy/Clone the alignment
```go
var clonealign align.Alignment
var err error
clonealign,err = align.Clone()
```
* Get the sequence having a specific name
```go
var sequence string
var err error
sequence,err = align.GetSequence("nameofsequence")
```
* Build a bootstrap replicate
```go
var bootstrap align.Alignment
bootstrap = align.BuildBootstrap()
```
* Randomly shuffle sequence order of alignment
```go
align.ShuffleSequences()
```
* Compute evolutionary ditance matrix (5 threads)
```go
import "github.com/fredericlemoine/goalign/distance"
//...
var model distance.DistModel
var distMatrix [][]float64
model = distance.Model("k2p", false)
distmatrix = distance.DistMatrix(align, nil, model, 5)
```

* Other functions

Other functions are described in the godoc.
