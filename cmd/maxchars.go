package cmd

import (
	"github.com/fredericlemoine/goalign/align"
	"github.com/fredericlemoine/goalign/io"
	"github.com/spf13/cobra"
)

// charCmd represents the char command
var maxCharCmd = &cobra.Command{
	Use:   "maxchar",
	Short: "Prints the character with the highest occcurence for each site of the alignment",
	Long: `Prints the character with the highest occcurence for each site of the alignment.

Ouput format: Tabulated with columns:
1) Site index (0...)
2) Character with maximum occurence
3) Number of occurence of this character

Example of usages:

goalign stats maxchar -i align.phylip -p
goalign stats maxchar -i align.fasta
`,
	RunE: func(cmd *cobra.Command, args []string) (err error) {
		var aligns *align.AlignChannel

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
		printMaxCharStats(al)

		return
	},
}

func init() {
	statsCmd.AddCommand(maxCharCmd)
}
