package cmd

import (
	"github.com/spf13/cobra"
)

var siteRate float64

// sitesCmd represents the sites command
var sitesCmd = &cobra.Command{
	Use:   "sites",
	Short: "Shuffles n alignment sites vertically",
	Long: `Shuffles n alignment sites vertically.

It may take Fasta or Phylip input alignment.

If the input alignment contains several alignments, will process all of them

The only option is to specify the rate of shuffled sites.
A rate of 0.5 will shuffle 50% of the sites of the alignment.

It takes n sites of the input alignment and reassign the 
characters to different sequences.

Example of usage:

goalign shuffle sites -i align.phylip -p -r 0.5
goalign shuffle sites -i align.fasta -r 0.5

`,
	Run: func(cmd *cobra.Command, args []string) {
		f := openWriteFile(shuffleOutput)
		for al := range rootaligns {
			al.ShuffleSites(siteRate)
			writeAlign(al, f)
		}
		f.Close()
	},
}

func init() {
	shuffleCmd.AddCommand(sitesCmd)
	sitesCmd.PersistentFlags().Float64VarP(&siteRate, "rate", "r", 0.5, "Rate of shuffled sites (>=0 and <=1)")
}
