package cmd

import (
	"fmt"
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
	Run: func(cmd *cobra.Command, args []string) {
		for al := range rootaligns {
			fmt.Println(al.Length())
		}
	},
}

func init() {
	statsCmd.AddCommand(lengthCmd)
}
