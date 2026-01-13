# Goalign: toolkit and api for alignment manipulation

## Commands

### subseq
This command takes a sub-alignment (slice) from the input alignment.

It takes an alignment and extracts sub-sequences from it, given
a start position (0-based inclusive) and a length.
If the length is after the end of the alignment, will stop at the 
end of the alignment.

(start,length) can be specified on the alignment coordinate system (default) 
or on a given sequence coordinate system without gaps (`--ref-seq`). It allows 
to extract sub alignment by just knowing coordinates on a given unaligned sequence.


If the length l == 0 , then the extracted sequences will be [start,alilength[
If the length l <  0 , then the extracted sequences will be [start,alilength-l[
If the length l <  0 and a reference sequence is given, the sub alignment will span [start,reflength-l[
of the ref sequence.

For example:
```
goalign subseq -p -i al.phy -s 9 -l 10
```

This will extract a sub-alignment going from 10th position of the alignemnt, with a length of 10.

The output format is the same than input format.

If `--ref-seq <name>` is specified, then the coordinates are specified on the given sequence
coordinate system without considering gaps.

For example:

If al.fa is:
```
>s1
--ACG--AT-GC
>s2
GGACGTTATCGC
```

Then:
```
goalign subseq -i al.fa -s 0 -l 4 --ref-seq s1
````

will output:

```
>s1
ACG--A
>s2
ACGTTA
```

`--ref-seq` is currently not not compatible with `--step`


Sliding window:
---------------

If `--step` is given and > 0, then Several sub-alignments will be produced,
and corresponding to all alignments in windows of sizes -l, and with starts:
[start, start+step, ..., end-length].

Example with an alignment al.phy of size 10 (0123456789)

```
goalign subseq -i al.phy -s 0 -l 5 --step 1
```
will produce alignments with positions:

```
01234
 12345
  23456
   34567
    45678
     56789
```

Warning: If output is stdout, it works only if input format is Phylip, because 
it is possible to split alignments afterwards (with `goalign divide` for example).

`--step` is currently not compatible with `--ref-seq`.

Several alignments:
------------------

If several alignments are present in the input file and the output is a file (not stdout or -) , then :

* First alignment, first subalignment: results will be placed in the given file
  (ex out.fasta)
* First alignment, other subalignments (sliding windows): results will be placed
  in file with the given name with `_sub<i>` suffix (ex: `out_sub1.fasta`, `out_sub2.fasta`, etc.)
* Other alignments, first subalignment: results will be placed in the given file
  with `_al<i>` suffix (ex `out_al1.fasta`, `out_al2.fasta`, etc.)
* Other alignments, other subalignments: results will be placed in the given file
  with `_al<i>` and `_sub<i>` suffixes (ex `out_al1_sub1.fasta`, `out_al1_sub2.fasta`, etc.)

In any case:
------------

If `--reverse` is given, all positions but the given subsequence is output. In case of a sliding window, it will correspond to sliding subsequence "deletion". In case of reference sequence coordinates, the subsequence is first computed, then all positions of the alignment but this subsequence are given in output.


#### Usage
```
Usage:
  goalign subseq [flags]
  
Flags:
  -l, --length int      Length of the sub alignment (default 10)
  -o, --output string   Alignment output file (default "stdout")
  -s, --start int       Start position
      --step int        Step: If > 0, then will generate several alignments, 
                        for each window of length l, with starts: 
                        [start,start+step, ..., end-l]*

Global Flags:
  -i, --align string   Alignment input file (default "stdin")
  -p, --phylip         Alignment is in phylip? False=Fasta
  --input-strict       Strict phylip input format (only used with -p)
  --output-strict      Strict phylip output format  (only used with -p)
```

#### Examples

* Generating a random alignment and taking a sub-alignment of length 10 from position 5:
```
goalign random -n 4 --seed 10 -l 100 | goalign subseq -l 10 -s 5
```

It should give the following alignment:
```
>Seq0000
ATTTGCCGTA
>Seq0001
GACCTAAGTA
>Seq0002
CGCGGGCGGA
>Seq0003
GAGGCTTTAT
```
