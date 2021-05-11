package cmd

import (
	"fmt"

	"github.com/evolbioinfo/goalign/align"
	"github.com/evolbioinfo/goalign/io"
	"github.com/spf13/cobra"
)

// nalignCmd represents the nalign command
var nalignCmd = &cobra.Command{
	Use:   "nalign",
	Short: "Prints the number of alignments in the input file",
	Long: `Prints the number of alignments in the input file

If the input file is in Fasta format, it should be 1
Otherwize, it may be > 1

Example:

goalign stats nalign -i align.ph -p

`,
	RunE: func(cmd *cobra.Command, args []string) (err error) {
		var aligns *align.AlignChannel

		if aligns, err = readalign(infile); err != nil {
			io.LogError(err)
			return
		}

		naligns := 0
		for range aligns.Achan {
			naligns++
		}
		fmt.Println(naligns)

		if aligns.Err != nil {
			err = aligns.Err
			io.LogError(err)
		}
		return
	},
}

func init() {
	statsCmd.AddCommand(nalignCmd)
}
