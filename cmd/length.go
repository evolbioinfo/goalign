package cmd

import (
	"fmt"

	"github.com/fredericlemoine/goalign/align"
	"github.com/fredericlemoine/goalign/io"
	"github.com/spf13/cobra"
)

// lengthCmd represents the length command
var lengthCmd = &cobra.Command{
	Use:   "length",
	Short: "Prints the length of sequences in the alignment",
	Long: `Prints the length of sequences in the alignment. 
May take a Phylip of Fasta input alignment.

If the input alignment contains several alignments, will take all of them

Example of usages:

goalign stats length -i align.phylip -p
goalign stats length -i align.fasta

`,
	RunE: func(cmd *cobra.Command, args []string) (err error) {
		var aligns align.AlignChannel

		if aligns, err = readalign(infile); err != nil {
			io.LogError(err)
			return
		}

		for al := range aligns.Achan {
			fmt.Println(al.Length())
		}

		if aligns.Err != nil {
			err = aligns.Err
			io.LogError(err)
		}
		return
	},
}

func init() {
	statsCmd.AddCommand(lengthCmd)
}
