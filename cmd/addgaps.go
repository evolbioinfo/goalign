package cmd

import (
	"github.com/spf13/cobra"
)

var gaprate float64
var gapnbseqs float64

// rogueCmd represents the rogue command
var addgapsCmd = &cobra.Command{
	Use:   "gaps",
	Short: "Adds random (uniformly) gaps",
	Long: `Adds random (uniformly) gaps on an input alignment.

Example:
goalign shuffle gaps -i align.fa -n 0.5 -r 0.5

`,
	Run: func(cmd *cobra.Command, args []string) {

		f := openWriteFile(shuffleOutput)
		nameFile := openWriteFile(rogueNameFile)
		for al := range rootaligns {
			al.AddGaps(gaprate, gapnbseqs)
			writeAlign(al, f)
		}
		nameFile.Close()
		f.Close()
	},
}

func init() {
	shuffleCmd.AddCommand(addgapsCmd)

	addgapsCmd.PersistentFlags().Float64VarP(&gapnbseqs, "prop-seq", "n", 0.5, "Proportion of the sequences to add gaps")
	addgapsCmd.PersistentFlags().Float64VarP(&gaprate, "gap-rate", "r", 0.5, "Proportion of gaps to add per sequences")
}
