package cmd

import (
	"github.com/spf13/cobra"
)

// sampleCmd represents the sample command
var sampleCmd = &cobra.Command{
	Use:   "sample",
	Short: "Samples sequences or sites from an input alignment",
	Long: `Samples sequences or sites from an input alignment. For example:

Randomly sampling 10 sequences from the alignment:
goalign sample seqs -i align.fa -n 10 > sample.fa

Randomly sampling 10 subsequences with length 20 from the input alignment:
goalign sample seqs -n 10 -i align.fa -l 20 -o subalign_

`,
}

func init() {
	RootCmd.AddCommand(sampleCmd)
}
