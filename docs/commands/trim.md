# Goalign: toolkit and api for alignment manipulation

## Commands

### trim
This command trims names of sequences or sequences themselves.

Two sub-commands:
* `goalign trim name`: trims sequence names to n characters. It will also output the correspondance between old names and new names into a map file as well as the new alignment.
* `goalign trim seq`: trims sequences from the left or from the right side, by n characters.

#### Usage
* General command:
```
Usage:
  goalign trim [command]

Available Commands:
  name        Trims names of sequences
  seq         Trims sequences of the alignment

Flags:
  -o, --out-align string   Trimed alignment output file (default "stdout")

Global Flags:
  -i, --align string   Alignment input file (default "stdin")
  -p, --phylip         Alignment is in phylip? False=Fasta
```

* `goalign trim name`: 
```
Usage:
  goalign trim name [flags]

Flags:
  -n, --nb-char int      Number of characters to keep in sequence names (default 1)
  -m, --out-map string   Mapping output file (default "none")

Global Flags:
  -i, --align string       Alignment input file (default "stdin")
  -o, --out-align string   Renamed alignment output file (default "stdout")
  -p, --phylip             Alignment is in phylip? False=Fasta
```


* `goalign trim seq`:
```
Usage:
  goalign trim seq [flags]

Flags:
  -s, --from-start    If true: trims n char from the left, otherwise from the right
  -n, --nb-char int   Number of characters to trim from sequences (default 1)

Global Flags:
  -i, --align string       Alignment input file (default "stdin")
  -o, --out-align string   Renamed alignment output file (default "stdout")
  -p, --phylip             Alignment is in phylip? False=Fasta
```

#### Examples
* Generating a random alignment and trimming sequence names :
```
goalign random -s 10 -n 4 -l 5 | goalign trim name -n 3 -m mapfile.txt
```

Should output:
```
>S01
GATTA
>S02
ATTTG
>S03
CCGTA
>S04
GGCCA
```
And print in `mapfile.txt`:
```
Seq0000	S01
Seq0001	S02
Seq0002	S03
Seq0003	S04
```

* Generating a random alignment and trimming 5 nt from the left:
```
goalign random -s 10 -n 4 -l 10 
goalign random -s 10 -n 4 -l 10 | goalign trim seq -n 5 -s 
```

Should output:
```
>Seq0000
GATTAATTTG
>Seq0001
CCGTAGGCCA
>Seq0002
GAATCTGAAG
>Seq0003
ATCGAACACT
```
Then:
```
>Seq0000
ATTTG
>Seq0001
GGCCA
>Seq0002
TGAAG
>Seq0003
ACACT
```
