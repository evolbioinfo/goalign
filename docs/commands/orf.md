# Goalign: toolkit and api for alignment manipulation

## Commands

### orf
Find the longest orf in all given sequences in forward strand.

If input sequences are not nucleotidic, then returns an error.

If input sequences are aligned (contain '-'), then they are unaligned first.

Output is in fasta format (format options such as -p and -x are ignored).

#### Usage
```
Usage:
  goalign orf [flags]

Flags:
  -h, --help            help for orf
  -o, --output string   ORF Output Fasta File (default "stdout")

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

