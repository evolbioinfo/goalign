# Goalign: toolkit and api for alignment manipulation

## Commands

### subset
Take a subset of sequences from the input alignment

It takes an alignment and a set of sequence names, and prints the alignments corresponding to sequence names

There are two ways of specifying sequence names:
1) Via `-f` option; file should be formated with one sequence name per line and or coma separated. If the file contains names that do not exist in the alignment, they won't be taken into account.
2) Via arguments on command line (ex: `goalign subset -i align.fa Seq0001 Seq0002`).

In both ways (`-f` or via command line arguments), if option `-e` is given, names are treated as regexps, and if option `--indices` is given they are treated as indices. For example, it is possible to keep only sequences whose name contains "human" with: `goalign subset -i align.fa -e "human"`, or with the case insensitive regexp : `goalign subset -i align.fa -e "(?i)human"`. In addition it is possible to keep only sequences with indices 2 and 349 with (3rd and 350th sequences) : `goalign subset --indices -i align.fa 2 349`

If several names/indices are given, it is considered a "OR", then to get all sequences whose name contains human OR mouse (case insensitive): `goalign subset -i align.fa -e "(?i)human" "(?i)mouse"`.

Finally, one can revert the matching with `-r` option. In that case, given sequences are removed instead.

subset may take unaligned sequences as input, in that case, --unaligned must be specified, and only fasta input format is accepted.

#### Usage
```
Usage:
  goalign subset [flags]
  
Flags:
  -h, --help               help for subset
  -f, --name-file string   File containing names of sequences to keep (default "stdin")
  -o, --output string      Alignment output file (default "stdout")
  -e, --regexp             If sequence names are given as regexp patterns
  -r, --revert             If true, will remove given sequences instead of keeping only them
      --unaligned          Considers input sequences as unaligned and fasta format (phylip, nexus,... options are ignored)

Global Flags:
  -i, --align string    Alignment input file (default "stdin")
      --auto-detect     Auto detects input format (overrides -p, -x and -u)
  -u, --clustal         Alignment is in clustal? default fasta
      --input-strict    Strict phylip input format (only used with -p)
  -x, --nexus           Alignment is in nexus? default fasta
      --no-block        Write Phylip sequences without space separated blocks (only used with -p)
      --one-line        Write Phylip sequences on 1 line (only used with -p)
      --output-strict   Strict phylip output format (only used with -p)
  -p, --phylip          Alignment is in phylip? default fasta
```

#### Examples

* Generating a random alignment and taking a subset of the sequences:
```
goalign random -n 40 --seed 10 -l 10 | goalign subset Seq0001 Seq0002
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
goalign random -n 40 --seed 10 -l 10 | goalign subset -f seqs
```

It should give the following alignment:
```
>Seq0001
CCGTAGGCCA
>Seq0002
GAATCTGAAG
```


* Generating a random alignment and keeping only sequences 0001 to 0009
```
goalign random -n 4000 --seed 10 -l 10 | goalign subset -e "Seq000[1-9]"
```
