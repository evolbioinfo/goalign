package cmd

import (
	"fmt"

	"github.com/evolbioinfo/goalign/align"
	"github.com/evolbioinfo/goalign/io"
	"github.com/spf13/cobra"
)

var compalign string

var identicalCmd = &cobra.Command{
	Use:   "identical",
	Short: "Assess whether the two alignments are identical",
	Long: `Assess whether the two alignments are identical.

`,
	RunE: func(cmd *cobra.Command, args []string) (err error) {
		var aligns, compaligns *align.AlignChannel

		if compalign == "none" {
			io.LogError(fmt.Errorf("no compared alignment has been given"))
			return
		}

		if aligns, err = readalign(infile); err != nil {
			io.LogError(err)
			return
		}

		if compaligns, err = readalign(compalign); err != nil {
			io.LogError(err)
			return
		}

		al := <-aligns.Achan
		comp := <-compaligns.Achan

		if aligns.Err != nil {
			err = aligns.Err
			io.LogError(err)
			return
		}
		if compaligns.Err != nil {
			err = compaligns.Err
			io.LogError(err)
			return
		}

		fmt.Println(al.Identical(comp))

		return
	},
}

func init() {
	RootCmd.AddCommand(identicalCmd)
	identicalCmd.PersistentFlags().StringVarP(&compalign, "compared", "c", "none", "Compared alignment file")
}
