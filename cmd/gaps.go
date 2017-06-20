package cmd

import (
	"math/rand"

	"github.com/spf13/cobra"
)

// mutateCmd represents the mutate command
var gapsCmd = &cobra.Command{
	Use:   "snvs",
	Short: "Adds substitutions uniformly in an input alignment",
	Long: `Adds substitutions uniformly in an input alignment.

      if rate <= 0 : does nothing
      if rate > 1  : then rate = 1
`,
	Run: func(cmd *cobra.Command, args []string) {
		rand.Seed(mutateSeed)
		f := openWriteFile(mutateOutput)
		for al := range rootaligns {
			al.Mutate(mutateRate)
			writeAlign(al, f)
		}
		f.Close()
	},
}

func init() {
	mutateCmd.AddCommand(gapsCmd)
}
