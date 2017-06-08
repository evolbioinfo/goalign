# Goalign: toolkit and api for alignment manipulation

## Commands

### subset
Take a subset of sequences from the input alignment

It takes an alignment and a set of sequence names, and prints the alignments corresponding to sequence names. 

There are two ways of specifying sequence names:
1) Via `-f` option; file should be formated with one sequence name per line and or coma separated. If the file contains names that do not exist in the alignment, they won't be taken into account.
2) Via arguments on command line (ex: `goalign subset -i align.fa Seq0001 Seq0002`).

In both ways (`-f` or via command line arguments), if option `-e` is given, names are treated as regexps. For example, it is possible to keep only sequences whose name contains "human" with: `goalign subset -i align.fa -e "human"`, or with the case insensitive regexp : `goalign subset -i align.fa -e "(?i)human"`. 

If several names are given, it is considered a "OR", then to get all sequences whose name contains human OR mouse (case insensitive): `goalign subset -i align.fa -e "(?i)human" "(?i)mouse"`.

Finally, one can revert the matching with `-r` option. In that case, given sequences are removed instead.

#### Usage
```
Usage:
  goalign subset [flags]
  
Flags:
	-f, --name-file string   File containing names of sequences to keep (default "stdin")
	-o, --output string      Alignment output file (default "stdout")
	-e, --regexp             If sequence names are given as regexp patterns
	-r, --revert             If true, will remove given sequences instead of keeping only them
		  
Global Flags:
	-i, --align string   Alignment input file (default "stdin")
	-p, --phylip         Alignment is in phylip? False=Fasta
	-t, --threads int    Number of threads (default 1)
```

#### Examples

* Generating a random alignment and taking a subset of the sequences:
```
goalign random -n 40 -s 10 -l 10 | goalign subset Seq0001 Seq0002
```

It should give the following alignment:
```
>Seq0001
CCGTAGGCCA
>Seq0002
GAATCTGAAG
```

* Generating a random alignment and taking a subset of the sequences using a text file:
```
echo "Seq0001" > seqs
echo "Seq0002" >> seqs
goalign random -n 40 -s 10 -l 10 | goalign subset -f seqs
```

It should give the following alignment:
```
>Seq0001
CCGTAGGCCA
>Seq0002
GAATCTGAAG
```
