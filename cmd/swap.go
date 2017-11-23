package cmd

import (
	"github.com/spf13/cobra"
)

var swapRate float64

// swapCmd represents the swap command
var swapCmd = &cobra.Command{
	Use:   "swap",
	Short: "Swap portion of sequences in the input alignment",
	Long: `Swap portion of sequences in the input alignment.
It may take Fasta or Phylip input alignment.

If the input alignment contains several alignments, will process all of them

The only option is to specify the rate of swap.
A rate of 0.5 will swap 25% of the sequences with 
other 25% of the sequences at a random position.

swap 50% of sequences:

s1 CCCCCCCCCCCCCC    s1 CCCCCCCCCCCCCC
s2 AAAAAAAAAAAAAA => s2 AAAAAATTTTTTTT
s3 GGGGGGGGGGGGGG    s3 GGGGGGGGGGGGGG
s4 TTTTTTTTTTTTTT    s4 TTTTTTAAAAAAAA

Example of usage:

goalign shuffle swap -i align.phylip -p -r 0.5
goalign shuffle swap -i align.fasta -r 0.5

`,
	Run: func(cmd *cobra.Command, args []string) {
		f := openWriteFile(shuffleOutput)
		for al := range rootaligns.Achan {
			al.Swap(swapRate)
			writeAlign(al, f)
		}
		f.Close()
	},
}

func init() {
	shuffleCmd.AddCommand(swapCmd)

	swapCmd.PersistentFlags().Float64VarP(&swapRate, "rate", "r", 0.5, "Rate of Swap sequences (>=0 and <=1)")
}
