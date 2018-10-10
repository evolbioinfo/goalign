package cmd

import (
	"fmt"

	"github.com/fredericlemoine/goalign/align"
	"github.com/fredericlemoine/goalign/io"
	"github.com/spf13/cobra"
)

// charCmd represents the char command
var alphabetCmd = &cobra.Command{
	Use:   "alphabet",
	Short: "Prints the alphabet detected for the input alignment",
	Long: `Prints  the alphabet detected for the input alignment.

`,
	RunE: func(cmd *cobra.Command, args []string) (err error) {
		var aligns align.AlignChannel

		if aligns, err = readalign(infile); err != nil {
			io.LogError(err)
			return
		}

		al, _ := <-aligns.Achan
		if aligns.Err != nil {
			err = aligns.Err
			io.LogError(err)
			return
		}
		fmt.Println(al.AlphabetStr())

		return
	},
}

func init() {
	statsCmd.AddCommand(alphabetCmd)
}
