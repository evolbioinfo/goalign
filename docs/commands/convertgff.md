# Goalign: toolkit and api for alignment manipulation

## Commands

### convertgff
This command converts the coordinates of an input GFF annotation from one sequence coordinate system in an alignment to another sequence coordinate system in the same alignment.

It reads an input alignment together with a GFF file that describes gene and CDS features. For each CDS feature, it uses the alignment to map the coordinates from the reference sequence to the destination sequence and writes a new GFF file with corresponding gene and CDS entries on the destination sequence.

The input GFF file is expected to contain gene lines and CDS lines. Each CDS line should be linked to its parent gene through the `Parent` attribute, and each gene line should define a `Name` attribute. The output file contains one gene line and one CDS line per converted feature.

The conversion is based on the alignment coordinates and the sequences specified by `--ref-seq` and `--dest-seq`. If a coordinate falls outside the alignment, or if a block length is invalid, the command exits with an error.

#### Usage
```
Usage:
  goalign convertgff [flags]

Flags:
      --dest-seq string   Destination sequence to which coordinates are converted (default "none")
      --gff string         GFF file with all coordinates to convert (default "none")
  -h, --help               help for convertgff
  -o, --output string      Output folder (default "stdout")
      --ref-seq string     Reference sequence on which coordinates are given (default "none")

Global Flags:
  -i, --align string          Alignment input file (default "stdin")
      --auto-detect           Auto detects input format (overrides -p, -x and -u)
  -u, --clustal               Alignment is in clustal? default fasta
      --input-strict          Strict phylip input format (only used with -p)
  -x, --nexus                 Alignment is in nexus? default fasta
      --no-block              Write Phylip sequences without space separated blocks (only used with -p)
      --one-line              Write Phylip sequences on 1 line (only used with -p)
      --output-strict         Strict phylip output format (only used with -p)
  -p, --phylip                Alignment is in phylip? default fasta
      --seed int              Random Seed: -1 = nano seconds since 1970/01/01 00:00:00 (default -1)
  -t, --threads int           Number of threads (default 1)
```

#### Examples
Convert coordinates from one reference sequence to another and write the result to a new GFF file:

```
goalign convertgff -i alignment.fasta --gff input.gff --ref-seq ref1 --dest-seq dest1 -o output.gff
```

This produces a GFF file where the gene and CDS coordinates are expressed on `dest1` instead of `ref1`.
