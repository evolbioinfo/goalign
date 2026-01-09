package cmd

import (
	"github.com/evolbioinfo/goalign/align"
	"github.com/evolbioinfo/goalign/io"
	"github.com/evolbioinfo/goalign/io/utils"
	"github.com/spf13/cobra"
)

// mutateCmd represents the mutate command
var addSNVsCmd = &cobra.Command{
	Use:   "snvs",
	Short: "Adds substitutions uniformly in an input alignment",
	Long: `Adds substitutions uniformly in an input alignment.

      if rate <= 0 : does nothing
      if rate > 1  : then rate = 1
`,
	RunE: func(cmd *cobra.Command, args []string) (err error) {
		var aligns *align.AlignChannel
		var f utils.StringWriterCloser

		if aligns, err = readalign(infile); err != nil {
			io.LogError(err)
			return
		}
		if f, err = utils.OpenWriteFile(mutateOutput); err != nil {
			io.LogError(err)
			return
		}
		defer utils.CloseWriteFile(f, mutateOutput)

		for al := range aligns.Achan {
			al.Mutate(mutateRate, globalRand)
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
	mutateCmd.AddCommand(addSNVsCmd)
}
