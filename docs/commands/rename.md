# Goalign: toolkit and api for alignment manipulation

## Commands

### rename
This command renames all sequences of the input alignment (fasta or phylip) using a map file. The map file  is tab separated, with the following fields:

1. Current name of sequence
2. New name of the sequence

If option `--revert` is given, it is the other way.

Note:
* If a sequence name does not appear in the map file, it will not be renamed;
* If a name that does not exist appears in the map file, it will not do anything;
* If input is in Phylip format with several alignments, all of them will be renamed.



#### Usage
```
Usage:
  goalign rename [flags]
  
  Flags:
    -m, --map-file string   Name Mapping infile (default "none")
	-o, --output string     renamed alignment output file (default "stdout")
	-r, --revert            Reverse orientation of mapfile
		
  Global Flags:
	-i, --align string   Alignment input file (default "stdin")
	-p, --phylip         Alignment is in phylip? False=Fasta
    --input-strict       Strict phylip input format (only used with -p)
    --output-strict      Strict phylip output format  (only used with -p)
```

#### Examples

* Generating a random alignment and renaming sequences:
```
echo -e "Seq0000\tNewSeq" > map
goalign random -n 3 -s 10 -l 10 | goalign rename -m map
```

It should give the following alignment:
```
>NewSeq
GATTAATTTG
>Seq0001
CCGTAGGCCA
>Seq0002
GAATCTGAAG
```
