# Goalign: toolkit and api for alignment manipulation

## Commands

### stats
This command prints different characteristics of the alignment.
Without any subcommand, it prints the following information:
1. Length of alignment;
2. Number of sequences;
3. Average number of alleles per site;
4. Number of variables sites (does ot take into account gaps or special characters);
5. Character frequencies.

Different sub-commands:
* `goalign stats alleles`: Prints the average number of alleles per site of the alignment;
* `goalign stats char`: Prints the character frequencies;
* `goalign stats length`: Prints alignment length;
* `goalign stats nalign`: Prints the number of alignments in the input file (Phylip);
* `goalign stats nseq`: Prints the number of sequences in the input alignment;
* `goalign stats taxa`: Lists taxa in the input alignment.

#### Usage
* General command:
```
Usage:
goalign stats [flags]
goalign stats [command]
  
Available Commands:
	alleles
	char
	length
	nalign
	nseq
	taxa
			  
Global Flags:
  -i, --align string   Alignment input file (default "stdin")
  -p, --phylip         Alignment is in phylip? False=Fasta
  --input-strict       Strict phylip input format (only used with -p)
```

#### Examples
* Generating a random (uniform) alignment and printing stats:
```
goalign random -l 20 -s10| goalign stats
```

Should give:
```
length 20
nseqs 10
avgalleles 3.7500
char nb freq
A 54 0.270000
C 34 0.170000
G 49 0.245000
T 63 0.315000
```
