package cmd

import (
	"github.com/fredericlemoine/goalign/io"
	"github.com/spf13/cobra"
	"math/rand"
	"time"
)

var sampleseqSeed int64
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
	Run: func(cmd *cobra.Command, args []string) {
		rand.Seed(sampleseqSeed)
		f := openWriteFile(sampleseqOutput)
		for al := range rootaligns.Achan {
			if sample, err := al.Sample(sampleseqNb); err != nil {
				io.ExitWithMessage(err)
			} else {
				writeAlign(sample, f)
			}
		}
		f.Close()
	},
}

func init() {
	sampleCmd.AddCommand(sampleseqCmd)
	sampleseqCmd.PersistentFlags().IntVarP(&sampleseqNb, "nb-seq", "n", 1, "Number of sequences to sample from the alignment")
	sampleseqCmd.PersistentFlags().Int64VarP(&sampleseqSeed, "seed", "s", time.Now().UTC().UnixNano(), "Initial Random Seed")
	sampleseqCmd.PersistentFlags().StringVarP(&sampleseqOutput, "output", "o", "stdout", "Sampled alignment output file")
}
