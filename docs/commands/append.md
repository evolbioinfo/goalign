# Goalign: toolkit and api for alignment manipulation

## Commands

### append
Append alignments to an input alignment by inserting new sequences.

This commands adds the sequences of a set of alignments to a reference alignement specified by `-i`.

If sequences do not have the same length than the reference alignment, then returns an error.

If format is phylip, it may contain several alignments in one file and then we can append all of them at once:

```
goalign append -i refalign.phy aligns.phy
```

If format is Fasta, several alignments may be given in the form:

```
goalign append -i align.fasta others*.fasta
```

#### Usage

General command

```
Usage:
  goalign append [flags]

Flags:
  -h, --help            help for append
  -o, --output string   Alignment output file (default "stdout")

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

* Append 3 alignments

```
cat > input.1 <<EOF
>s1
ACGACGACGACC
>2
ATCTT-TTTTTC
>3
ATCTT-TTTTTT
EOF

cat > input.2 <<EOF
>s4
ACGACGACGACC
>5
ATCTT-TTTTTC
>6
ATCTT-TTTTTT
EOF

cat > input.3 <<EOF
>s7
ACGACGACGACC
>8
ATCTT-TTTTTC
>9
ATCTT-TTTTTT
EOF

goalign append -i input.1 input.2 input.3 

>s1
ACGACGACGACC
>2
ATCTT-TTTTTC
>3
ATCTT-TTTTTT
>s4
ACGACGACGACC
>5
ATCTT-TTTTTC
>6
ATCTT-TTTTTT
>s7
ACGACGACGACC
>8
ATCTT-TTTTTC
>9
ATCTT-TTTTTT

```
