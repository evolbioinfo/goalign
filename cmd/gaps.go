package cmd

import (
	"os"

	"github.com/fredericlemoine/goalign/align"
	"github.com/fredericlemoine/goalign/io"
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
	RunE: func(cmd *cobra.Command, args []string) (err error) {
		var aligns align.AlignChannel
		var f *os.File

		if aligns, err = readalign(infile); err != nil {
			io.LogError(err)
			return
		}
		if f, err = openWriteFile(mutateOutput); err != nil {
			io.LogError(err)
			return
		}
		defer closeWriteFile(f, mutateOutput)

		for al := range aligns.Achan {
			al.Mutate(mutateRate)
			writeAlign(al, f)
		}

		if aligns.Err != nil {
			err = aligns.Err
			io.LogError(err)
		}
		return
	},
}

func init() {
	mutateCmd.AddCommand(gapsCmd)
}
