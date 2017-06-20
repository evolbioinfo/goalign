package cmd

import (
	"math/rand"

	"github.com/spf13/cobra"
)

var gapnbseqs float64

// rogueCmd represents the rogue command
var addgapsCmd = &cobra.Command{
	Use:   "gaps",
	Short: "Adds gaps uniformly in an input alignment",
	Long: `Adds gaps uniformly in an input alignment.

Example:
goalign mutate gaps -i align.fa -n 0.5 -r 0.5
`,
	Run: func(cmd *cobra.Command, args []string) {
		rand.Seed(mutateSeed)
		f := openWriteFile(mutateOutput)
		for al := range rootaligns {
			al.AddGaps(mutateRate, gapnbseqs)
			writeAlign(al, f)
		}
		f.Close()
	},
}

func init() {
	mutateCmd.AddCommand(addgapsCmd)
	addgapsCmd.PersistentFlags().Float64VarP(&gapnbseqs, "prop-seq", "n", 0.5, "Proportion of the sequences in which to add gaps")
}
