package cmd

import (
	"github.com/evolbioinfo/goalign/align"
	"github.com/evolbioinfo/goalign/io"
	"github.com/evolbioinfo/goalign/io/utils"
	"github.com/spf13/cobra"
)

var nnbseqs float64

// rogueCmd represents the rogue command
var addAmbigCmd = &cobra.Command{
	Use:   "ambig",
	Short: "Adds ambiguities uniformly to an input alignment",
	Long: `Adds ambiguities uniformly to an input alignment.

Example:
goalign mutate ambig -i align.fa -n 0.5 -r 0.5
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
			al.AddAmbiguities(mutateRate, gapnbseqs)
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
	mutateCmd.AddCommand(addAmbigCmd)
	addAmbigCmd.PersistentFlags().Float64VarP(&nnbseqs, "prop-seq", "n", 0.5, "Proportion of the sequences in which to add ambiguities")
}
