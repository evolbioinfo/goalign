# Goalign: toolkit and api for alignment manipulation

## Commands

### compress
This command removes identical patterns/sites from an input alignment

#### Usage
```
Usage:
  goalign compress [flags]

Flags:
  -o, --output string       Compressed output alignment file (default "stdout")
      --weight-out string   Pattern weight output file (default "none")

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

```
cat > input <<EOF
   4   20
1 GGGGGGGGGGGGGGGGGGGG
2 TTTTTTTTTTTTTTTTTTTT
3 GGGGGGGGGGCCCCTTTTTT
4 AAAAAAAAAAAAAAAAAAAA
EOF
```

goalign compress -i input -p --weight-out wres

should output:
```
   4   3
1 GGG
2 TTT
3 CGT
4 AAA
```

and produce a weight file like:
```
4
10
6
```
