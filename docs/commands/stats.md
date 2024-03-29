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

If `--per-sequences` option is given, then prints statistics for all sequences individually:
1. sequence: Sequence ID;
2. gaps: Total number of gaps;
3. gapsstart: Number of consecutive gaps  at the beginning;
4. gapsend: Number of consecutive gaps at the end;
5. gapsuniques: Number of gaps unique to that sequence;

	If --count-profile is given, then the output will be : unique\tnew\tboth, with:
	- 5a unique: # gaps that are unique in each sequence in the alignment
	- 5b new: # gaps that are new in each sequence compared to the profile
	- 5c both: # gaps that are unique in each sequence in the alignment and that are new compared the profile

6. gapsopenning: Number of streches of gaps;
7. mutuniques: Number of unique mutations;

	If --count-profile is given, then the output will be : unique\tnew\tboth, with:
	- 7a unique: # mutations that are unique in each sequence in the alignment
	- 7b new: # mutations that are new in each sequence compared to the profile
	- 7c both: # mutations that are unique in each sequence in the alignment and that are new compared the profile

7. (bis) mutref:  Number mutations compared to a reference sequence (only if `--ref-sequence`is given)
8. length: Lenght of the unaligned sequence;
9. A	C	G	T...: Number of occurence of each character.

Note that `--count-profile`takes a tab separated file such as given by the command `goalign stats char --per-sites`:

```
site  A C G T
0 nA  nC  nG  nT
1...
...
n
```

Different sub-commands:
* `goalign stats alleles`: Prints the average number of alleles per site of the alignment;
* `goalign stats alphabet`: Prints the alphabet of the alignemnts (aminoacids, nucleotides, unknown);
* `goalign stats char`: Prints the character number of occurences. If `--per-sequences` is given, then prints the number of occurences of each characters for each seqences. If `--per-sites` is given, then prints the number of occurences of each characters for each sites. Is is possible to give `--only` option, to count the number of occurences of a single character.
* `goalign stats gaps`: Prints the number of gaps in each sequences (and possibly the number of gaps from start, and from end); By default, it prints, for each alignment sequence the number of gaps. Following options are exclusive, and given in order of priority: If `--from-start` is specified, then counts only gaps at sequence starts; If `--from-end` is specified, then counts only gaps at sequence ends; If `--unique` is specified, then counts only gaps that are unique in their alignmebnnt columnIf` --openning` is specified, then counts only gap openning (streches of gaps are counted once); Otherwise, counts total number of gaps on each sequence. If `--profile` is given in addition to `--unique`, then the output will be : `unique\tnew\tboth`, with:

  - unique: # gaps that are unique in their column, for each sequence of the alignment
  - new: # gaps that are new in each sequence compared to the profile
  - both: # gaps that are unique in each sequence in the alignment and that are new compared the profile.

* `goalign stats length`: Prints alignment length;
* `goalign stats maxchar`: Prints max occurence char for each alignment site. It can ignore gaps or Ns with `--ignore-gaps` and `--ignore-n` (N/n for nucleotides and X/x if amino-aciods) (except if only Ns or Gaps at the position);
* `goalign stats mutations`: Prints, for each sequence, the number of mutations on each alignment sequence, compared to a reference sequence or unique compared to all other sequences. It does not take into account '-' and 'N' as unique mutations, and does not take into account '-' and 'N' as mutations compared to a reference sequence; 	If `--unique` is specified, then counts only mutations (characters) that are unique in their column for the given sequence.	If `--ref-sequence` is specified, it will try to extract a sequence having that name from the alignment. If none exist, it will try to open a fasta file with the given name to take the first sequence as a reference. If a character is ambigous (IUPAC notation) in an nucleotide sequence, then it is counted as a mutation only if it is incompatible with the reference character. If `--profile` is given in addition to `--unique`, then the output will be : `unique\tnew\tboth`, with:

  - unique: # mutations that are unique in their column, for each sequence of the alignment
  - new: # mutations that are new in each sequence compared to the profile
  - both: # mutations that are unique in each sequence in the alignment and that are new compared the profile.
	It does not take into account 'N' as mutations compared to a reference sequence.
  
* `goalign stats mutations list` : Print mutation list of each alignment sequence compared to the given reference sequence. Options:
	- --ref-sequence: it will try to extract the given sequence from the alignment. If none exist, 
	it will try to open a fasta file with the given name to take the first sequence as a reference. If a character is ambigous 
	(IUPAC notation) in a nucleotide sequence, then it is counted as a mutation only if it is incompatible with the reference character.
	- --aa : takes reference sequence codon by codon to list mutations in the aligned sequences. In case of an insertion or a deletion in the target sequence: if length%3!=0 (without gaps): it may be a frameshift, indicated by a '/'. It is better to use this option rather than translating the alignment and then listing mutations in aa, because the insertions/deletions may not be appropriately listed if the gap is inside a reference codon for example.


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
  -i, --align string          Alignment input file (default "stdin")
      --auto-detect           Auto detects input format (overrides -p, -x and -u)
  -u, --clustal               Alignment is in clustal? default fasta
      --ignore-identical int  Ignore duplicated sequences that have the same name and same sequences
      --input-strict          Strict phylip input format (only used with -p)
  -x, --nexus                 Alignment is in nexus? default fasta
      --no-block              Write Phylip sequences without space separated blocks (only used with -p)
      --one-line              Write Phylip sequences on 1 line (only used with -p)
      --output-strict         Strict phylip output format (only used with -p)
  -p, --phylip                Alignment is in phylip? default fasta
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
