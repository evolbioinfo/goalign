# Goalign: toolkit and api for alignment manipulation

## Commands

### subsites
This command takes a subset of the sites of the given alignment from the input alignment.

User provides a set of positions to extract from the input alignment, either in the command line or in an input file (one position per line, 0-based coordinates).

Site position (0-based, inclusive) can be specified on the alignment coordinate system (default) 
or on a given sequence coordinate system without gaps (`--ref-seq`). It allows 
to extract positions by just knowing coordinates on a given unaligned reference sequence.

For example:
```
goalign subsites -p -i al.phy 1 2 3
```

or

```
goalign subsites -p -i al.phy --sitefile positions.txt
```


The output format is the same than input format.

If `--ref-seq <name>` is specified, then the coordinates are specified on the given sequence
coordinate system without considering gaps.

For example:

If al.fa is:
```
>s1
--ACG--AT-GC
>s2
GGACGTTATCGC
```

Then:
```
goalign subsites -i al.fa --ref-seq s1 0 1 2 3
````

will output:

```
>s1
ACGA
>s2
ACGA
```

If `--informative` is given, only informative sites (parsimony definition) are selected. Informative sites are the positions that contain at least 2 different characters that occur at least twice each. This option has priority over the index based site selection above, `--ref-seq` is ignored and `--reverse` is still taken into account.

In any case:
------------

If `--reverse` is given, all positions but the given positions are output. In case of reference sequence coordinates, the positions are first extracted based on reference coordinates, then all positions of the alignment but these are given in output.


#### Usage
```
Usage:
  goalign subsites [flags]

Flags:
  -h, --help              help for subsites
  --informative           Selects (~parsimony) informative sites
  -o, --output string     Alignment output file (default "stdout")
      --ref-seq string    Reference sequence on which coordinates are given (default "none")
  -r, --reverse           Take all but the given sites
      --sitefile string   File with positions of sites to select (one perline) (default "none")

Global Flags:
  -i, --align string       Alignment input file (default "stdin")
      --auto-detect        Auto detects input format (overrides -p, -x and -u)
  -u, --clustal            Alignment is in clustal? default fasta
      --ignore-identical   Ignore duplicated sequences that have the same name and same sequences
      --input-strict       Strict phylip input format (only used with -p)
  -x, --nexus              Alignment is in nexus? default fasta
      --no-block           Write Phylip sequences without space separated blocks (only used with -p)
      --one-line           Write Phylip sequences on 1 line (only used with -p)
      --output-strict      Strict phylip output format (only used with -p)
  -p, --phylip             Alignment is in phylip? default fasta
```


