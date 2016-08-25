package cmd

import (
	"github.com/fredericlemoine/goalign/io/phylip"
	"github.com/spf13/cobra"
)

// phylipCmd represents the phylip command
var phylipCmd = &cobra.Command{
	Use:   "phylip",
	Short: "Reformats an input alignment into Phylip",
	Long: `Reformats an alignment into Phylip. 
It may take a Phylip of Fasta input alignment.

Example of usage:

goalign reformat phylip -i align.phylip -p
goalign reformat phylip -i align.fasta

`,
	Run: func(cmd *cobra.Command, args []string) {
		reformatOutputString = phylip.WriteAlignment(rootalign)
	},
}

func init() {
	reformatCmd.AddCommand(phylipCmd)
}
