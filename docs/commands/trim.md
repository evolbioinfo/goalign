# Goalign: toolkit and api for alignment manipulation

## Commands

### trim
This command trims names of sequences or sequences themselves.

Two sub-commands:
* `goalign trim name`: trims sequence names to n characters. It will also output the correspondance between old names and new names into a map file as well as the new alignment. If `-a` is given, then generates sequence names automatically. If `--unaligned` is given, sequences are considered unaligned.
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
  -a, --auto             Automatically generates sequence identifiers (priority over --nb-cchar)
  -h, --help             help for name
  -n, --nb-char int      Number of characters to keep in sequence names (default 1)
  -m, --out-map string   Mapping output file (default "none")
      --unaligned        Considers sequences as unaligned and format fasta (phylip, nexus,... options are ignored)

Global Flags:
  -i, --align string           Alignment input file (default "stdin")
      --auto-detect            Auto detects input format (overrides -p, -x and -u)
  -u, --clustal                Alignment is in clustal? default fasta
      --ignore-identical int   Ignore duplicated sequences that have the same name and potentially have same sequences, 0 : Does not ignore anything, 1: Ignore sequences having the same name (keep the first one whatever their sequence), 2: Ignore sequences having the same name and the same sequence
      --input-strict           Strict phylip input format (only used with -p)
  -x, --nexus                  Alignment is in nexus? default fasta
      --no-block               Write Phylip sequences without space separated blocks (only used with -p)
      --one-line               Write Phylip sequences on 1 line (only used with -p)
  -o, --out-align string       Renamed alignment output file (default "stdout")
      --output-strict          Strict phylip output format (only used with -p)
  -p, --phylip                 Alignment is in phylip? default fasta
```


* `goalign trim seq`:
```
Usage:
  goalign trim seq [flags]

Flags:
  -s, --from-start    If true: trims n char from the start, else from the end
  -h, --help          help for seq
  -n, --nb-char int   Number of characters to trim from sequences (default 1)

Global Flags:
  -i, --align string           Alignment input file (default "stdin")
      --auto-detect            Auto detects input format (overrides -p, -x and -u)
  -u, --clustal                Alignment is in clustal? default fasta
      --ignore-identical int   Ignore duplicated sequences that have the same name and potentially have same sequences, 0 : Does not ignore anything, 1: Ignore sequences having the same name (keep the first one whatever their sequence), 2: Ignore sequences having the same name and the same sequence
      --input-strict           Strict phylip input format (only used with -p)
  -x, --nexus                  Alignment is in nexus? default fasta
      --no-block               Write Phylip sequences without space separated blocks (only used with -p)
      --one-line               Write Phylip sequences on 1 line (only used with -p)
  -o, --out-align string       Renamed alignment output file (default "stdout")
      --output-strict          Strict phylip output format (only used with -p)
  -p, --phylip                 Alignment is in phylip? default fasta
      --seed int               Random Seed: -1 = nano seconds since 1970/01/01 00:00:00 (default -1)
  -t, --threads int            Number of threads (default 1)
```

#### Examples
* Generating a random alignment and trimming sequence names :
```
goalign random --seed 10 -n 4 -l 5 | goalign trim name -n 3 -m mapfile.txt
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
goalign random --seed 10 -n 4 -l 10 
goalign random --seed 10 -n 4 -l 10 | goalign trim seq -n 5 -s 
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
