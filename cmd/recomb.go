package cmd

import (
	"github.com/spf13/cobra"
)

var recombRate float64

// recombCmd represents the recomb command
var recombCmd = &cobra.Command{
	Use:   "recomb",
	Short: "Recombines sequences in the input alignment",
	Long: `Recombines sequences in the input alignment.
It may take Fasta or Phylip input alignment.

The only option is to specify the rate of recombination.
A rate of 0.5 will recombine 25% of the sequences with 
other 25% of the sequences at a random position.

Example of usage:

goalign shuffle recombine -i align.phylip -p -r 0.5
goalign shuffle recombine -i align.fasta -r 0.5

`,
	Run: func(cmd *cobra.Command, args []string) {
		rootalign.Recombine(recombRate)
	},
}

func init() {
	shuffleCmd.AddCommand(recombCmd)

	recombCmd.PersistentFlags().Float64VarP(&recombRate, "rate", "r", 0.5, "Rate of recombined sequences (>=0 and <=1)")
}
