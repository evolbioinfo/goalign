# Goalign: toolkit and api for alignment manipulation

## Commands

### replace
This command replaces characters in sequences of an input alignment. The `--regexp (-e)` option considers the string to be replaced as a regular expression.
Unless `--unaligned`is specified, the replacement should not change sequence lengths, otherwise it returns an error.

If `--posfile` is given, then `--old` and `--new` are not considered. Instead, characters at sites+sequences specified in the input file
are replaced in the alignement. The format of the input posfile is tabulated with columns:

- 0: sequence name
- 1: site index
- 2: new character

#### Usage
```
Usage:
  goalign replace [flags]

Flags:
  -h, --help            help for replace
  -n, --new string      New string that will replace old string in sequences (default "none")
  -s, --old string      String to replace in the sequences (default "none")
  -o, --output string   Output alignment file (default "stdout")
  -e, --regexp          Considers Replace alignment using regexp
  -f, --posfile string  File containing sites to replace by give characters in given sequences (deactivates --old & --new) (default "none")
      --unaligned       Considers input sequences as unaligned and fasta format (phylip, nexus,... options are ignored)

Global Flags:
  -i, --align string    Alignment input file (default "stdin")
      --auto-detect     Auto detects input format (overrides -p, -x and -u)
  -u, --clustal         Alignment is in clustal? default fasta
      --input-strict    Strict phylip input format (only used with -p)
  -x, --nexus           Alignment is in nexus? default fasta
      --no-block        Write Phylip sequences without space separated blocks (only used with -p)
      --one-line        Write Phylip sequences on 1 line (only used with -p)
      --output-strict   Strict phylip output format (only used with -p)
  -p, --phylip          Alignment is in phylip? default fasta
```

#### Examples

* Replacing all `` with `---` (regex):

Input alignment:
```
   10   20
Seq0000  GATTAATTTG CCGTAGGCCA
Seq0001  GAATCTGAAG ATCGAACACT
Seq0002  TTAAGTTTTC ACTTCTAATG
Seq0003  GAGAGGACTA GTTCATACTT
Seq0004  TTTAAACACT TTTACATCGA
Seq0005  TGTCGGACCT AAGTATTGAG
Seq0006  TACAACGGTG TATTCCAGCG
Seq0007  GTGGAGAGGT CTATTTTTCC
Seq0008  GGTTGAAGGA CTCTAGAGCT
Seq0009  GTAAAGGGTA TGGCCATGTG
```

```
goalign replace -s 'GA.' -n '---' -e -p -i input -o result
```

It should give the following alignment:
```
   10   20
Seq0000  ---TAATTTG CCGTAGGCCA
Seq0001  ---TCT---- --C---CACT
Seq0002  TTAAGTTTTC ACTTCTAATG
Seq0003  ---AG---TA GTTCATACTT
Seq0004  TTTAAACACT TTTACATCGA
Seq0005  TGTCG---CT AAGTATT---
Seq0006  TACAACGGTG TATTCCAGCG
Seq0007  GTG---AGGT CTATTTTTCC
Seq0008  GGTT---G-- -TCTA---CT
Seq0009  GTAAAGGGTA TGGCCATGTG
```

### replace stops

Replace STOP codons in input nt sequences by NNN, in the given phase (except the last codon).

If --phase is given (>=0), then starts at the given offset (default=0)

```
Usage:
  goalign replace stops [flags]

Flags:
      --genetic-code string   Genetic Code: standard, mitoi (invertebrate mitochondrial) or mitov (vertebrate mitochondrial) (default "standard")
  -h, --help                  help for stops
      --phase int             Phase in which replace STOP codons

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
  -o, --output string          Output alignment file (default "stdout")
      --output-strict          Strict phylip output format (only used with -p)
  -p, --phylip                 Alignment is in phylip? default fasta
  -k, --stockholm              Alignment is in stockholm? default fasta
      --unaligned              Considers input sequences as unaligned and fasta format (phylip, nexus,... options are ignored)
```
