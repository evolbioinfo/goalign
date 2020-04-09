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
6. Aphabet

Different sub-commands:
* `goalign stats alleles`: Prints the average number of alleles per site of the alignment;
* `goalign stats alphabet`: Prints the alphabet of the alignemnts (aminoacids, nucleotides, unknown);
* `goalign stats char`: Prints the character occurences;
* `goalign stats gaps`: Prints the number of gaps in each sequences (and possibly the number of gaps from start, and from end);
* `goalign stats length`: Prints alignment length;
* `goalign stats maxchar`: Prints max occurence char for each alignment site;
* `goalign stats mutations`: Prints the number of mutations on each alignment sequence, compared to a reference sequence;
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
  alleles     Prints the average number of alleles per sites of the alignment
  alphabet    Prints the alphabet detected for the input alignment
  char        Prints frequence of different characters (aa/nt) of the alignment
  gaps        Print gap stats on each alignment sequence
  length      Prints the length of sequences in the alignment
  maxchar     Prints the character with the highest occcurence for each site of the alignment
  mutations   Print mutations stats on each alignment sequence compared to a reference sequence
  nalign      Prints the number of alignments in the input file
  nseq        Prints the number of sequences in the alignment
  taxa        Prints index (position) and name of taxa of the alignment file
			  
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

#### Examples
* Generating a random (uniform) alignment and printing stats:
```
goalign random -l 20 --seed 10| goalign stats
```

Should give:
```
length	20
nseqs	10
avgalleles	3.7500
variable sites	20
char	nb	freq
A	54	0.270000
C	34	0.170000
G	49	0.245000
T	63	0.315000
alphabet	nucleotide
```
