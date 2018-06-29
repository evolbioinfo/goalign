# Goalign: toolkit and api for alignment manipulation

## Commands

### rename
This command renames all sequences of the input alignment (fasta or phylip) by 2 ways:

* Using a map file. The map file  is tab separated, with the following fields:

  1. Current name of sequence
  2. New name of the sequence

  If option `--revert` is given, it is the other way.

  Note:
    * If a sequence name does not appear in the map file, it will not be renamed;
    * If a name that does not exist appears in the map file, it will not do anything;
    * If input is in Phylip format with several alignments, all of them will be renamed.

* Using a regexp (`-e|--regexp`) and a replacement string (`-b|--replace`):
   will replace matching strings in sequence names by string given by `-b`, ex: `goalign rename -i align.fasta --regexp 'Seq(\d+)' --replace 'Newname$1' -m map.txt`
  this will replace all matches of 'Seq(\d+)' with 'NewName$1', with `$1` being the matched string inside ().


#### Usage
```
Usage:
  goalign rename [flags]

Flags:
  -h, --help              help for rename
  -m, --map-file string   Name Mapping infile (default "none")
  -o, --output string     renamed alignment output file (default "stdout")
  -e, --regexp string     rename alignment using given regexp (default "none")
  -b, --replace string    replaces regexp matching strings by this string (default "none")
  -r, --revert            Reverse orientation of mapfile

Global Flags:
  -i, --align string    Alignment input file (default "stdin")
      --auto-detect     Auto detects input format (overrides -p and -x)
      --input-strict    Strict phylip input format (only used with -p)
  -x, --nexus           Alignment is in nexus? default fasta
      --output-strict   Strict phylip output format (only used with -p)
  -p, --phylip          Alignment is in phylip? default fasta
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

* Generating a random alignment and renaming sequences with a regex:
```
goalign random -n 3 -s 10 -l 10 | goalign rename --regexp 'Seq' --replace 'NewSeq' -m map.txt
```

It should give the following alignment:
```
>NewSeq0000
GATTAATTTG
>NewSeq0001
CCGTAGGCCA
>NewSeq0002
GAATCTGAAG
```
