# Goalign: toolkit and api for alignment manipulation

## Commands

### mask
Mask a part of the input alignment (replace by N|X or a specified character)

It takes an alignment and replaces some characters by "N" or "X" by default.

By default, it masks positions defined by start (0-based inclusive)
and length with N (nucleotide alignment) or X (amino-acid alignment) by default.
If the length is after the end of the alignment, will stop at the 
end of the alignment.

For example:
```
goalign mask -p -i al.phy -s 9 -l 10
```

This will replace 10 positions with N|X from the 10th position.

If `--ref-seq` is specified, then coordinates are considered on the given reference sequence
without considering gaps. In that case, if the range of masked sites incorporates gaps in
the reference sequence, these gaps will also be masked in the output alignment.

If `--unique` is specified, `goalign mask --unique` will replace characters that
are unique (defined by `--at-most` option) in their column (except GAPS) with N or X. 
In that case, if `--ref-seq` option is given, then a unique character is masked if:
    1) It is different from the given reference sequence
    2) Or the reference is a GAP
In this case, `--length` and `--start` are ignored.

Option `--replace` defines the replacement character. If `--replace` is "" (default), then, 
masked characters are replaced by "N" or "X" depending on the alphabet. 
Orherwise:
  1) if `--replace` is AMBIG: just like ""
  2) if `--replace` is MAJ: Masked characters are replaced by the most frequent character of the column 
     (without considering the reference sequence when `--ref-seq` is given)
  3) if `--replace` is GAP: Replacing character is a GAP

The output format is the same than input format.

#### Usage
```
Usage:
  goalign mask [flags]

Flags:
  -h, --help            help for mask
  -l, --length int      Length of the sub alignment (default 10)
  -o, --output string   Alignment output file (default "stdout")
  -s, --start int       Start position (0-based inclusive)
      --unique          If given, then masks characters that are unique in their columns (start and length are ignored)

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

align.phy :
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
goalign mask -i align.phy -p -s 0 -l 2
```

Should print:
```
   10   20
Seq0000  NNTTAATTTG CCGTAGGCCA
Seq0001  NNATCTGAAG ATCGAACACT
Seq0002  NNAAGTTTTC ACTTCTAATG
Seq0003  NNGAGGACTA GTTCATACTT
Seq0004  NNTAAACACT TTTACATCGA
Seq0005  NNTCGGACCT AAGTATTGAG
Seq0006  NNCAACGGTG TATTCCAGCG
Seq0007  NNGGAGAGGT CTATTTTTCC
Seq0008  NNTTGAAGGA CTCTAGAGCT
Seq0009  NNAAAGGGTA TGGCCATGTG
```
