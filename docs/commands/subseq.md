# Goalign: toolkit and api for alignment manipulation

## Commands

### subseq
This command takes a sub-alignment (slice) from the input alignment.

It takes an alignment and extracts sub-sequences from it, given a start position (0-based inclusive) and a length. If the length is after the end of the alignment, will stop at the end of the alignment.

#### Usage
```
Usage:
  goalign subseq [flags]
  
Flags:
	-l, --length int      Length of the sub alignment (default 10)
	-o, --output string   Alignment output file (default "stdout")
	-s, --start int       Start position
		
Global Flags:
	-i, --align string   Alignment input file (default "stdin")
	-p, --phylip         Alignment is in phylip? False=Fasta
    --input-strict       Strict phylip input format (only used with -p)
    --output-strict      Strict phylip output format  (only used with -p)
```

#### Examples

* Generating a random alignment and taking a sub-alignment of length 10 from position 5:
```
goalign random -n 4 -s 10 -l 100 | goalign subseq -l 10 -s 5
```

It should give the following alignment:
```
>Seq0000
ATTTGCCGTA
>Seq0001
GACCTAAGTA
>Seq0002
CGCGGGCGGA
>Seq0003
GAGGCTTTAT
```
