# Goalign: toolkit and api for alignment manipulation
## Github repository
[Goalign github repository](https://github.com/evolbioinfo/goalign).
## Introduction
Goalign is a set of command line tools to manipulate multiple alignments. It is implemented in [Go](https://golang.org/) language.

The goal is to handle multiple alignments in different input and output formats (Fasta, Phylip, Clustal and Nexus) through several basic commands. Each command may print result (usually an alignment) in the standard output, and thus can be piped to the standard input of the next Goalign command.

## Installation
### Binaries
You can download already compiled binaries for the latest release in the [release](https://github.com/evolbioinfo/goalign/releases) section.
Binaries are available for MacOS, Linux, and Windows (32 and 64 bits).

Once downloaded, you can just run the executable without any other downloads.

### From sources
In order to compile goalign, you must first [download](https://golang.org/dl/) and [install](https://golang.org/doc/install) Go on your system.

Then you just have to type :
```bash
go get github.com/evolbioinfo/goalign/
```
This will download Goalign sources from github, and all its dependencies.

You can then build it with:
```bash
cd $GOPATH/src/github.com/evolbioinfo/goalign/
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
* `--alphabet`: Used to specify which alphabet must be used to parse the alignment. It can be `auto` (default), `aa`, or `nt`. By default, the alphabet is deduced from the content of the input file. In the case of nexus format, when `--alphabet auto` is specified, the alphabet specified in the nexus file is used. Otherwise, this option overrides the nexus file alphabet.

Output is written to stdout by default, but can be generally written to files with the `-o` option. If the given output file has a `.gz` or `.xz` extension, the output is compressed accordingly.

Command                                                     | Subcommand |        Description
------------------------------------------------------------|------------|-----------------------------------------------------------------------
[addid](commands/addid.md) ([api](api/addid.md))            |            | Adds a string to each sequence identifier of the input alignment
[append](commands/append.md) ([api](api/append.md))         |            | Concatenates several alignments by adding new alignments as new sequences of the first alignment
[build](commands/build.md) ([api](api/build.md))            |            | Command to build output files : bootstrap for example
--                                                          | distboot   | Builds bootstrap distances matrices from input alignment (nt only)
--                                                          | seqboot    | Builds bootstrap alignments from input alignment
[clean](commands/clean.md) ([api](api/clean.md))            |            | Removes gap sites/sequences
--                                                          | sites      | Removes sequences with gaps
--                                                          | seqs       | Removes sites with gaps
[codonalign](commands/codonalign.md) ([api](api/codonalign.md))|         | Adds gaps in nt sequences, according to its corresponding protein alignment
[compress](commands/compress.md) ([api](api/compress.md))   |            | Removes identical patterns/sites from an input alignment
[compute](commands/compute.md) ([api](api/compute.md))      |            | Different computations (distances, entropy, etc.)
--                                                          | distance   | Computes distance matrix from inpu alignment
--                                                          | entropy    | Computes entropy of sites of a given alignment
--                                                          | pssm       | Computes and prints a Position specific scoring matrix
[concat](commands/concat.md) ([api](api/concat.md))         |            | Concatenates a set of alignment
[consensus](commands/consensus.md) ([api](api/consensus.md))|            | Computes a basic majority consensus sequence
[extract](commands/extract.md)                              |            | Extracts sub-sequences from an input alignment
[completion](commands/completion.md)                        |            | Generates auto-completion commands for bash or zsh
[dedup](commands/dedup.md) ([api](api/dedup.md))            |            | Deduplicate/Remove identical sequences 
[diff](commands/diff.md) ([api](api/diff.md))               |            | Compares all sequences of an alignment to the first one, and counts differences
[divide](commands/divide.md) ([api](api/divide.md))         |            | Divide an input alignment in several output files
[draw](commands/draw.md) ([api](api/draw.md))               |            | Draws an input alignment
--                                                          | biojs      | Displays an input alignment in an html file using biojs
--                                                          | png        | Displays an input alignment in a png file
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
[rename](commands/rename.md) ([api](api/rename.md))         |            | Rename sequences of the input alignment (using a map file, with a regexp, or just clean names)
[replace](commands/replace.md) ([api](api/replace.md))      |            | Replace characters in sequences of input alignment
--                                                          | stops      | Replace stop codons in input nt alignment
[revcomp](commands/revcomp.md) ([api](api/revcomp.md))      |            | Reverse complements an input alignment
[sample](commands/sample.md) ([api](api/sample.md))         |            | Samples sequences or sites from an input alignment
--                                                          | seqs       | Samples a subset of sequences from the input alignment
--                                                          | sites      | Takes a random subalignment
--                                                          | rarefy     | Takes a sample taking into accounts weights
[shuffle](commands/shuffle.md) ([api](api/shuffle.md))      |            | A set of commands to shuffle an alignment
--                                                          | recomb     | Recombines sequences in the input alignment (copy/paste)
--                                                          | rogue      | Simulates rogue taxa
--                                                          | seqs       | Shuffles sequence order in alignment
--                                                          | sites      | Shuffles n alignment sites vertically
--                                                          | swap       | Swaps portion of sequences in the input alignment (cut/paste)
[split](commands/split.md) ([api](api/split.md))            |            | Split an input alignment according to partitions defined in an partition file
[sort](commands/sort.md) ([api](api/sort.md))               |            | Sorts the alignment by sequence name
[stats](commands/stats.md) ([api](api/stats.md))            |            | Prints different characteristics of the alignment
--                                                          | alleles    | Prints the average number of alleles per sites of the alignment
--                                                          | alphabet   | Prints  the alphabet detected for the alignment
--                                                          | char       | Prints frequence of different characters (aa/nt) of the alignment
--                                                          | gaps       | Prints statistics about gaps for each sequence of the alignment
--                                                          | length     | Prints the length of sequences in the alignment
--                                                          | mutations  | Prints, for each sequence, the number of mutations compared to a reference sequence
--                                                          | maxchar    | Prints max occurence char for each alignment site
--                                                          | nalign     | Prints the number of alignments in the input file (phylip)
--                                                          | nseq       | Prints the number of sequences in the alignment
--                                                          | taxa       | Prints index (position) and name of taxa of the alignment file
[subseq](commands/subseq.md) ([api](api/subseq.md))         |            | Take a sub-alignment from the input alignment
[subset](commands/subset.md) ([api](api/subset.md))         |            | Take a subset of sequences from the input alignment
[subsites](commands/subsites.md) (api)                      |            | Take a subset of the sites from the input alignment
[sw](commands/sw.md) ([api](api/sw.md))                     |            | Aligns 2 sequences using Smith&Waterman algorithm
[tolower](commands/tolower.md) ([api](api/tolower.md))      |            | Replace upper case characters by lower case characters
[toupper](commands/toupper.md) ([api](api/toupper.md))      |            | Replace lower case characters by upper case characters
[translate](commands/translate.md) ([api](api/translate.md))|            | Translates an input sequence into Amino-Acids
[transpose](commands/transpose.md) ([api](api/transpose.md))|            | Transposes an input alignment (sequences<=>sites)
[trim](commands/trim.md) ([api](api/trim.md))               |            | This command trims names of sequences or sequences themselves
--                                                          | name       | Trims names of sequences
--                                                          | seq        | Trims sequences of the input alignment
[unalign](commands/unalign.md) ([api](api/unalign.md))      |            | Unaligns input alignment
[version](commands/version.md)                              |            | Prints the current version of goalign
