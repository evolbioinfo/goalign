package cmd

import (
	"bufio"
	goio "io"

	"github.com/evolbioinfo/goalign/align"
	"github.com/evolbioinfo/goalign/io"
	"github.com/evolbioinfo/goalign/io/fasta"
	"github.com/evolbioinfo/goalign/io/utils"
	"github.com/spf13/cobra"
)

var codonAlignOutput string
var nucleotideFasta string

// codonAlignCmd
var codonAlignCmd = &cobra.Command{
	Use:   "codonalign",
	Short: "Aligns a given nt fasta file using a corresponding aa alignment",
	Long: `Aligns a given nt fasta file using a corresponding aa alignment.

If the input alignment is not amino acid, then returns an error.
If the given fasta file is not nucleotides then returns an error.

Warning: It does not check that the amino acid sequence is a good 
translation of the nucleotide sequence, but just add gaps to the
nucleotide sequence where needed. 

Once gaps are added, if the nucleotide alignment length does not match 
the protein alignment length * 3, returns an error.
`,
	RunE: func(cmd *cobra.Command, args []string) (err error) {
		var aligns *align.AlignChannel
		var f utils.StringWriterCloser
		var ntseqsf *bufio.Reader
		var toclose goio.Closer
		var ntseqs align.SeqBag
		var codonAl align.Alignment

		// Read input aa alignment
		if aligns, err = readalign(infile); err != nil {
			io.LogError(err)
			return
		}

		// Read input fasta nt sequences
		if toclose, ntseqsf, err = utils.GetReader(nucleotideFasta); err != nil {
			io.LogError(err)
			return
		}
		defer toclose.Close()

		if ntseqs, err = fasta.NewParser(ntseqsf).ParseUnalign(); err != nil {
			io.LogError(err)
			return
		}

		// Open output file
		if f, err = utils.OpenWriteFile(codonAlignOutput); err != nil {
			io.LogError(err)
			return
		}
		defer utils.CloseWriteFile(f, codonAlignOutput)

		for al := range aligns.Achan {
			if codonAl, err = al.CodonAlign(ntseqs); err != nil {
				io.LogError(err)
				return
			}
			writeAlign(codonAl, f)
		}

		if aligns.Err != nil {
			err = aligns.Err
			io.LogError(err)
		}
		return
	},
}

func init() {
	RootCmd.AddCommand(codonAlignCmd)
	codonAlignCmd.PersistentFlags().StringVarP(&codonAlignOutput, "output", "o", "stdout", "Output codon aligned file")
	codonAlignCmd.PersistentFlags().StringVarP(&nucleotideFasta, "fasta", "f", "stdin", "Input nucleotide Fasta file to be codon aligned")
}
