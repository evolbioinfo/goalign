# Goalign: toolkit and api for alignment manipulation
## Github repository
[Goalign github repository](https://github.com/fredericlemoine/goalign).
## Introduction
Goalign is a set of command line tools to manipulate multiple alignments. It is implemented in [Go](https://golang.org/) language.

The goal is to handle multiple alignments in different input and output formats (Fasta, Phylip, Clustal and Nexus) through several basic commands. Each command may print result (usually an alignment) in the standard output, and thus can be piped to the standard input of the next Goalign command.

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
* `-u`: input is in clustan format (default fasta), lower priority than `-p` and `-x`. Output format will also be clustal in this case;
* `--input-strict`: if `-p` is also given, then input is considered phylip strict, i.e:
    * sequence names are maximum 10 character long. goalign removes spaces in sequence names;
	* sequence starts at position 11 (just after sequence name).
* `--output-strict`: if `-p` is also given, then output alignments are written in strict phylip format, i.e:
    * sequence names are maximum 10 character long, otherwise they are truncated;
* `--no-block`: if `-p`is also given, then output alignments are written in phylip, without 10 character block separation.
* `--one-line`: if `-p`is also given, then output alignments are written inphylip, on one single line.
* `--auto-detect` (overrides `-p`, `-u` and `-x`): It will test input formats in the following order:
    1. Fasta
    2. Nexus
	3. Clustal
    4. Phylip
    If none of these formats is recognized, then will exit with an error. Please also note that in `--auto-detect` mode, phylip format is considered as not strict.

Command                                                     | Subcommand |        Description
------------------------------------------------------------|------------|-----------------------------------------------------------------------
[addid](commands/addid.md) ([api](api/addid.md))            |            | Adds a string to each sequence identifier of the input alignment
[build](commands/build.md) ([api](api/build.md))            |            | Command to build output files : bootstrap for example
--                                                          | distboot   | Builds bootstrap distances matrices from input alignment (nt only)
--                                                          | seqboot    | Builds bootstrap alignments from input alignment
[clean](commands/clean.md) ([api](api/clean.md))            |            | Removes gap sites/sequences
--                                                          | sites      | Removes sequences with gaps
--                                                          | seqs       | Removes sites with gaps
[codonalign](commands/codonalign.md) ([api](api/codonalign.md))|         | Adds gaps in nt sequences, according to its corresponding protein alignment
[compute](commands/compute.md) ([api](api/compute.md))      |            | Different computations (distances, entropy, etc.)
--                                                          | distance   | Computes distance matrix from inpu alignment
--                                                          | entropy    | Computes entropy of sites of a given alignment
--                                                          | pssm       | Computes and prints a Position specific scoring matrix
[concat](commands/concat.md) ([api](api/concat.md))         |            | Concatenates a set of alignment
[completion](commands/completion.md)                        |            | Generates auto-completion commands for bash or zsh
[dedup](commands/dedup.md) ([api](api/dedup.md))            |            | Deduplicate/Remove identical sequences 
[divide](commands/divide.md) ([api](api/divide.md))         |            | Divide an input alignment in several output files
[draw](commands/draw.md) ([api](api/draw.md))               |            | Draws an input alignment
--                                                          | biojs      | Displays an input alignment in an html file using biojs
[identical](commands/identical.md) ([api](api/identical.md))|            | Tells whether two alignments are identical
[mask](commands/mask.md) ([api](api/mask.md))               |            | Mask (with N or X) positions of input alignment
[mutate](commands/mutate.md) ([api](api/mutate.md))         |            | Adds substitutions (~sequencing errors), or gaps, uniformly in an input alignment
--                                                          | gaps       | Adds gaps uniformly in an input alignment
--                                                          | snvs       | Adds substitutions uniformly in an input alignment
[orf](commands/orf.md) ([api](api/orf.md))                  |            | Find the longest orf in all given sequences in forward strand
[phase](commands/phase.md) ([api](api/phase.md))            |            | Find best Starts by aligning to translated ref sequences and set them as new start positions
[phasent](commands/phasent.md) ([api](api/phase.md))        |            | Find best Starts by aligning to ref sequences and set them as new start positions
[random](commands/random.md) ([api](api/random.md))         |            | Generate random sequences
[reformat](commands/reformat.md) ([api](api/reformat.md))   |            | Reformats input alignment into phylip of fasta format
--                                                          | clustal    | Reformats an input alignment into Clustal
--                                                          | fasta      | Reformats an input alignment into Fasta
--                                                          | nexus      | Reformats an input alignment into nexus
--                                                          | paml       | Reformats an input alignment into PAML input format
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
--                                                          | alphabet   | Prints  the alphabet detected for the alignment
--                                                          | char       | Prints frequence of different characters (aa/nt) of the alignment
--                                                          | length     | Prints the length of sequences in the alignment
--                                                          | maxchar    | Prints max occurence char for each alignment site
--                                                          | nalign     | Prints the number of alignments in the input file (phylip)
--                                                          | nseq       | Prints the number of sequences in the alignment
--                                                          | taxa       | Prints index (position) and name of taxa of the alignment file
[subseq](commands/subseq.md) ([api](api/subseq.md))         |            | Take a sub-alignment from the input alignment
[subset](commands/subset.md) ([api](api/subset.md))         |            | Take a subset of sequences from the input alignment
[sw](commands/sw.md) ([api](api/sw.md))                     |            | Aligns 2 sequences using Smith&Waterman algorithm
[translate](commands/translate.md) ([api](api/translate.md))|            | Translates an input sequence into Amino-Acids
[trim](commands/trim.md) ([api](api/trim.md))               |            | This command trims names of sequences or sequences themselves
--                                                          | name       | Trims names of sequences
--                                                          | seq        | Trims sequences of the input alignment
[unalign](commands/unalign.md) ([api](api/unalign.md))      |            | Unaligns input alignment
[version](commands/version.md)                              |            | Prints the current version of goalign
