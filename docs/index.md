# Goalign: toolkit and api for alignment manipulation
## Github repository
[Goalign github repository](https://github.com/fredericlemoine/goalign).
## Introduction
Goalign is a set of command line tools to manipulate multiple alignments. It is implemented in [Go](https://golang.org/) language.

The goal is to handle multiple alignments in different input and output formats (Fasta, Phylip and Nexus) through several basic commands. Each command may print result (usually an alignment) in the standard output, and thus can be piped to the standard input of the next Goalign command.

## Installation
### Binaries
You can download already compiled binaries for the latest release in the [release](https://github.com/fredericlemoine/goalign/releases) section.
Binaries are available for MacOS, Linux, and Windows (32 and 64 bits).

Once downloaded, you can just run the executable without any other downloads.

### From sources
In order to compile goalign, you must first [download](https://golang.org/dl/) and [install](https://golang.org/doc/install) Go on your system.

Then you just have to type :
```bash
go get github.com/fredericlemoine/goalign/
```
This will download Goalign sources from github, and all its dependencies.

You can then build it with:
```bash
cd $GOPATH/src/github.com/fredericlemoine/goalign/
make
```
The `goalign` executable should be located in the `$GOPATH/bin` folder.

## Commands

Here is the list of all commands, with the link to the full description, and a link to a snippet that does it in GO.
Almost all commands can have the following arguments:

* `-p`: input is in phylip format (default fasta). Output format will also be phylip in this case;
* `-x`: input is in nexus format (default fasta), lower priority than `-p`. Output format will also be nexus in this case;
* `--input-strict`: if `-p` is also given, then input is considered phylip strict, i.e:
    * sequence names are maximum 10 character long. goalign removes spaces in sequence names;
	* sequence starts at position 11 (just after sequence name).
* `--output-strict`: if `-p` is also given, then output alignments are written in strict phylip format, i.e:
    * sequence names are maximum 10 character long, otherwise they are truncated.

Command                                                     | Subcommand |        Description
------------------------------------------------------------|------------|-----------------------------------------------------------------------
[addid](commands/addid.md) ([api](api/addid.md))            |            | Adds a string to each sequence identifier of the input alignment
[build](commands/build.md) ([api](api/build.md))            |            | Command to build output files : bootstrap for example
--                                                          | distboot   | Builds bootstrap distances matrices from input alignment (nt only)
--                                                          | seqboot    | Builds bootstrap alignments from input alignment
[clean](commands/clean.md) ([api](api/clean.md))            |            | Removes gap sites/sequences
--                                                          | sites      | Removes sequences with gaps
--                                                          | seqs       | Removes sites with gaps
[compute](commands/compute.md) ([api](api/compute.md))      |            | Different computations (distances, entropy, etc.)
--                                                          | distance   | Computes distance matrix from inpu alignment
--                                                          | entropy    | Computes entropy of sites of a given alignment
--                                                          | pssm       | Computes and prints a Position specific scoring matrix
[concat](commands/concat.md) ([api](api/concat.md))         |            | Concatenates a set of alignment
[divide](commands/divide.md) ([api](api/divide.md))         |            | Divide an input alignment in several output files
[mutate](commands/mutate.md) ([api](api/mutate.md))         |            | Adds substitutions (~sequencing errors), or gaps, uniformly in an input alignment
--                                                          | gaps       | Adds gaps uniformly in an input alignment
--                                                          | snvs       | Adds substitutions uniformly in an input alignment
[random](commands/random.md) ([api](api/random.md))         |            | Generate random sequences
[reformat](commands/reformat.md) ([api](api/reformat.md))   |            | Reformats input alignment into phylip of fasta format
--                                                          | fasta      | Reformats an input alignment into Fasta
--                                                          | nexus      | Reformats an input alignment into nexus
--                                                          | phylip     | Reformats an input alignment into Phylip
--                                                          | tnt        | Reformats an input alignment into TNT input file
[rename](commands/rename.md) ([api](api/rename.md))         |            | Rename sequences of the input alignment, given a map file
[sample](commands/sample.md) ([api](api/sample.md))         |            | Samples sequences or sites from an input alignment
--                                                          | seqs       | Samples a subset of sequences from the input alignment
--                                                          | sites      | Take a random subalignment
--                                                          | rarefy     | Take a sample taking into accounts weights
[shuffle](commands/shuffle.md) ([api](api/shuffle.md))      |            | A set of commands to shuffle an alignment
--                                                          | recomb     | Recombines sequences in the input alignment (copy/paste)
--                                                          | rogue      | Simulates rogue taxa
--                                                          | seqs       | Shuffles sequence order in alignment
--                                                          | sites      | Shuffles n alignment sites vertically
--                                                          | swap       | Swaps portion of sequences in the input alignment (cut/paste)
[sort](commands/sort.md) ([api](api/sort.md))               |            | Sorts the alignment by sequence name
[stats](commands/stats.md) ([api](api/stats.md))            |            | Prints different characteristics of the alignment
--                                                          | alleles    | Prints the average number of alleles per sites of the alignment
--                                                          | char       | Prints frequence of different characters (aa/nt) of the alignment
--                                                          | length     | Prints the length of sequences in the alignment
--                                                          | nalign     | Prints the number of alignments in the input file (phylip)
--                                                          | nseq       | Prints the number of sequences in the alignment
--                                                          | taxa       | Prints index (position) and name of taxa of the alignment file
[subseq](commands/subseq.md) ([api](api/subseq.md))         |            | Take a sub-alignment from the input alignment
[subset](commands/subset.md) ([api](api/subset.md))         |            | Take a subset of sequences from the input alignment
[trim](commands/trim.md) ([api](api/trim.md))               |            | This command trims names of sequences or sequences themselves
--                                                          | name       | Trims names of sequences
--                                                          | seq        | Trims sequences of the input alignment
[unalign](commands/unalign.md) ([api](api/unalign.md))      |            | Unaligns input alignment
[version](commands/version.md)                              |            | Prints the current version of goalign
