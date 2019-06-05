# Goalign: toolkit and api for alignment manipulation

## Commands

### dedup
This command deduplicates sequences. To do so it removes identical sequences.

If -l is specified, then identical sequences are printed in the given file
with the following format:

```
seq1,seq2
seq3,seq4
```

This means that seq1 is identical to seq2 and seq3 is identical to seq4.

#### Usage
```
Usage:
  goalign dedup [flags]

Flags:
  -h, --help            help for dedup
  -l, --log string      Deduplicated output log file (default "none")
  -o, --output string   Deduplicated output alignment file (default "stdout")

Global Flags:
  -i, --align string    Alignment input file (default "stdin")
      --auto-detect     Auto detects input format (overrides -p, -x and -u)
  -u, --clustal         Alignment is in clustal? default fasta
      --input-strict    Strict phylip input format (only used with -p)
  -x, --nexus           Alignment is in nexus? default fasta
      --output-strict   Strict phylip output format (only used with -p)
  -p, --phylip          Alignment is in phylip? default fasta
```

#### Examples

```
cat > input.phy <<EOF
   6   6
3  GGGGGG
1  AAAAAA
2  CCCCCC
4  GGGGGG
4  GGGGGG
EOF

goalign dedup -i input.phy --phylip
```

should output:
```
   3   6
3  GGGGGG
1  AAAAAA
2  CCCCCC
```
