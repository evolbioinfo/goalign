package cmd

import (
	"github.com/spf13/cobra"

	"github.com/fredericlemoine/goalign/io"
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
	Run: func(cmd *cobra.Command, args []string) {
		aligns := readalign(infile)
		al, _ := <-aligns.Achan
		if aligns.Err != nil {
			io.ExitWithMessage(aligns.Err)
		}
		printMaxCharStats(al)
	},
}

func init() {
	statsCmd.AddCommand(maxCharCmd)
}
