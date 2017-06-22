# Goalign: toolkit and api for alignment manipulation

## Commands

### reformat
This command reformats an input alignment (fasta by default or phylip with `-p`) in different formats depending on the sub-command:
1. `goalign reformat fasta`: reformats input alignment in fasta;
2. `goalign reformat nexus`: reformats input alignment in nexus;
3. `goalign reformat phylip`: reformats input alignment in phylip;
4. `goalign reformat tnt`: reformats input alignment in TNT input format.


#### Usage
* general command:
```
Usage:
  goalign reformat [command]

Available Commands:
  fasta       Reformats an input alignment into Fasta
  nexus       Reformats an input alignment into nexus
  phylip      Reformats an input alignment into Phylip
  tnt         Reformats an input alignment into input data for TNT

Flags:
  -o, --output string   Reformated alignment output file (default "stdout")
  -s, --strict          If it is strict phylip format

Global Flags:
  -i, --align string   Alignment input file (default "stdin")
  -p, --phylip         Alignment is in phylip? False=Fasta
  --input-strict       Strict phylip input format (only used with -p)
  --output-strict      Strict phylip output format  (only used with -p)
```

#### Examples

* Generating a random fasta alignment with 5 sequences and 5 nucleotides and reformating it in phylip:
```
goalign random -n 5 -l 5 -s 10 | goalign reformat phylip
```

Should give the following alignment:
```
  5   5
Seq0000  GATTA
Seq0001  ATTTG
Seq0002  CCGTA
Seq0003  GGCCA
Seq0004  GAATC
```
* Generating a random fasta alignment with 5 sequences and 5 nucleotides and reformating it in nexus:
```
goalign random -n 5 -l 5 -s 10 | goalign reformat nexus
```

Should give the following alignment:
```
#NEXUS
begin data;
dimensions ntax=5 nchar=5;
format datatype=nucleotide;
matrix
Seq0000 GATTA
Seq0001 ATTTG
Seq0002 CCGTA
Seq0003 GGCCA
Seq0004 GAATC
;
end;
```
