# Goalign: toolkit and api for alignment manipulation

## Commands

### toupper
This command replaces lower case characters by upper case characters.

#### Usage
```
Usage:
  goalign toupper [flags]

Flags:
  -h, --help            help for tolower
  -o, --output string   Output upper case alignment file (default "stdout")
      --unaligned       Considers sequences as unaligned and format fasta (phylip, nexus,... options are ignored)

Global Flags:
  -i, --align string           Alignment input file (default "stdin")
      --alphabet string        Alignment/Sequences alphabet: auto (default), aa, or nt (default "auto")
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

```
echo ">Seq0000
gattaatttg
>Seq0001
ccgtaggcca
>Seq0002
gaatctgaag
>Seq0003
atcgaacact
>Seq0004
ttaagttttc" | ./goalign toupper
```

It should give the following alignment:
```
>Seq0000
GATTAATTTG
>Seq0001
CCGTAGGCCA
>Seq0002
GAATCTGAAG
>Seq0003
ATCGAACACT
>Seq0004
TTAAGTTTTC
```
