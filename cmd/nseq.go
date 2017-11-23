package cmd

import (
	"fmt"

	"github.com/spf13/cobra"
)

// nseqCmd represents the nseq command
var nseqCmd = &cobra.Command{
	Use:   "nseq",
	Short: "Prints the number of sequences in the alignment",
	Long: `Prints the number of sequences in the alignment. 
May take a Phylip of Fasta input alignment.

If the input alignment contains several alignments, will process all of them

Example of usages:

goalign stats nseq -i align.phylip -p
goalign stats nseq -i align.fasta
`,
	Run: func(cmd *cobra.Command, args []string) {
		for al := range rootaligns.Achan {
			fmt.Println(al.NbSequences())
		}
	},
}

func init() {
	statsCmd.AddCommand(nseqCmd)
}
