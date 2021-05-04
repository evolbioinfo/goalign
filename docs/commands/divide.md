# Goalign: toolkit and api for alignment manipulation

## Commands

### divide
The default behavior is to take an input alignment file containing potentially several alignments (e.g. with Phylip format ), and utput one alignment per output file.

If the alignment is in fasta format : will create 1 file. Otherwise, will create one output file per input alignment.

If the option `--nb-sequences <n>` is given, then will print n sequences per output file.

`-o` is the prefix of output files

Ex: if `-o div`, it will create files div_000.ph...div_n.ph

Output files will be in the same format as input files, or in fasta if `-f` is given.

Example:

gotree divide -i align.ph -p -o out


#### Usage
```
Usage:
  goalign divide [flags]

Flags:
  -h, --help               help for divide
      --nb-sequences int   Number of sequences per output file (<=0 : all sequences, >0: each alignment will be additionnaly split in several alignments)
  -f, --out-fasta          Forces output files to be in fasta format
  -o, --output string      Divided alignment output files prefix (default "prefix")
      --unaligned          Considers sequences as unaligned and format fasta (phylip, nexus,... options are ignored)

Global Flags:
  -i, --align string           Alignment input file (default "stdin")
      --auto-detect            Auto detects input format (overrides -p, -x and -u)
  -u, --clustal                Alignment is in clustal? default fasta
      --ignore-identical int   Ignore duplicated sequences that have the same name and potentially have same sequences, 0 : Does not ignore anything, 1: Ignore sequences having the same name (keep the first one whatever their sequence), 2: Ignore sequences having the same name and the same sequence
      --input-strict           Strict phylip input format (only used with -p)
  -x, --nexus                  Alignment is in nexus? default fasta
      --no-block               Write Phylip sequences without space separated blocks (only used with -p)
      --one-line               Write Phylip sequences on 1 line (only used with -p)
      --output-strict          Strict phylip output format (only used with -p)
  -p, --phylip                 Alignment is in phylip? default fasta
```

#### Examples

Input file, input.phylip: 
```
   11   10
Seq0000  GATTAATTTG
Seq0001  CCGTAGGCCA
Seq0002  GAATCTGAAG
Seq0003  ATCGAACACT
Seq0004  TTAAGTTTTC
Seq0005  ACTTCTAATG
Seq0006  GAGAGGACTA
Seq0007  GTTCATACTT
Seq0008  TTTAAACACT
Seq0009  TTTACATCGA
Seq0010  TGTCGGACCT
   3   10
Seq0000  GATTAATTTG
Seq0001  CCGTAGGCCA
Seq0002  GAATCTGAAG
```

* `goalign divide -i input -p -o divprefix -f --nb-sequences 2` will generate the following files:

* `divprefix_000.fa`
```
>Seq0000
GATTAATTTG
>Seq0001
CCGTAGGCCA
```
* `divprefix_001.fa`
```
>Seq0002
GAATCTGAAG
>Seq0003
ATCGAACACT
```
* `divprefix_002.fa`
```
>Seq0004
TTAAGTTTTC
>Seq0005
ACTTCTAATG
```
* `divprefix_003.fa`
```
>Seq0006
GAGAGGACTA
>Seq0007
GTTCATACTT
```
* `divprefix_004.fa`
```
>Seq0008
TTTAAACACT
>Seq0009
TTTACATCGA
```
* `divprefix_005.fa`
```
>Seq0010
TGTCGGACCT
```
* `divprefix_006.fa`
```
>Seq0000
GATTAATTTG
>Seq0001
CCGTAGGCCA
```
* `divprefix_007.fa`
```
>Seq0002
GAATCTGAAG
```
