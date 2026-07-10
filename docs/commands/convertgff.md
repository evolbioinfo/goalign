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
goalign convertgff -i input.fasta --gff input.gff --ref-seq seq1 --dest-seq seq2 -o output.gff
```

Example alignment file `input.fasta`:

```fasta
>seq1
AC--ACGTACGT
>seq2
ACGT----ACGT
>seq3
ACGTACGT--GT
```

Example GFF file `input.gff`:

```gff
seq1	Genbank	gene	1	4	.	+	.	ID=gene-g1;Name=Gene1;gene=N;gene_biotype=protein_coding
seq1	Genbank	CDS	1	4	.	+	.	ID=cds-g1;Name=CDS1;gene=N;gene_biotype=protein_coding;Parent=gene-g1
seq1	Genbank	gene	5	8	.	+	.	ID=gene-g2;Name=Gene2;gene=P;gene_biotype=protein_coding
seq1	Genbank	CDS	5	8	.	+	.	ID=cds-g2;Name=CDS2;gene=P;gene_biotype=protein_coding;Parent=gene-g2
seq1	Genbank	gene	1	6	.	+	.	ID=gene-g3;Name=Gene3;gene=Q;gene_biotype=protein_coding
seq1	Genbank	CDS	1	6	.	+	.	ID=cds-g3;Name=CDS3;gene=Q;gene_biotype=protein_coding;Parent=gene-g3
```

This produces a GFF file where the gene and CDS coordinates are expressed on `seq2` instead of `seq1`, for example:

```gff
seq2	Genbank	gene	5	6	.	+	.	ID=gene-Gene2;Name=Gene2;gbkey=Gene;gene=Gene2;gene_biotype=protein_coding
seq2	Genbank	CDS	5	6	.	+	.	ID=cds-Gene2;Name=Gene2;gbkey=CDS;Parent=gene-Gene2
seq2	Genbank	gene	1	4	.	+	.	ID=gene-Gene3;Name=Gene3;gbkey=Gene;gene=Gene3;gene_biotype=protein_coding
seq2	Genbank	CDS	1	4	.	+	.	ID=cds-Gene3;Name=Gene3;gbkey=CDS;Parent=gene-Gene3
seq2	Genbank	gene	1	4	.	+	.	ID=gene-Gene1;Name=Gene1;gbkey=Gene;gene=Gene1;gene_biotype=protein_coding
seq2	Genbank	CDS	1	4	.	+	.	ID=cds-Gene1;Name=Gene1;gbkey=CDS;Parent=gene-Gene1
```
