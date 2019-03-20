package cmd

import (
	"os"

	"github.com/fredericlemoine/goalign/align"
	"github.com/fredericlemoine/goalign/io"
	"github.com/spf13/cobra"
)

var sampleseqOutput string
var sampleseqNb int

// sampleCmd represents the sample command
var sampleseqCmd = &cobra.Command{
	Use:   "seqs",
	Short: "Samples a subset of sequences from the input alignment",
	Long: `Samples a subset of sequences from the input alignment.

May take a Fasta or Phylip alignment in input.

If the input alignment contains several alignments, will process all of them

As output, writes an alignment containing a sample of the sequences

`,
	RunE: func(cmd *cobra.Command, args []string) (err error) {
		var aligns *align.AlignChannel
		var f *os.File
		var sample align.Alignment

		if aligns, err = readalign(infile); err != nil {
			io.LogError(err)
			return
		}
		if f, err = openWriteFile(sampleseqOutput); err != nil {
			io.LogError(err)
			return
		}
		defer closeWriteFile(f, sampleseqOutput)

		for al := range aligns.Achan {
			if sample, err = al.Sample(sampleseqNb); err != nil {
				io.LogError(err)
				return
			} else {
				writeAlign(sample, f)
			}
		}

		if aligns.Err != nil {
			err = aligns.Err
			io.LogError(err)
		}
		return
	},
}

func init() {
	sampleCmd.AddCommand(sampleseqCmd)
	sampleseqCmd.PersistentFlags().IntVarP(&sampleseqNb, "nb-seq", "n", 1, "Number of sequences to sample from the alignment")
	sampleseqCmd.PersistentFlags().StringVarP(&sampleseqOutput, "output", "o", "stdout", "Sampled alignment output file")
}
