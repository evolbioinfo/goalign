package cmd

import (
	"github.com/fredericlemoine/goalign/io/fasta"
	"github.com/spf13/cobra"
)

// fastaCmd represents the fasta command
var fastaCmd = &cobra.Command{
	Use:   "fasta",
	Short: "Reformats an input alignment into Fasta",
	Long: `Reformats an alignment into Fasta. 
It may take a Phylip of Fasta input alignment.

Example of usage:

goalign reformat fasta -i align.phylip -p
goalign reformat fasta -i align.fasta

`,
	Run: func(cmd *cobra.Command, args []string) {
		reformatOutputString = fasta.WriteAlignment(rootalign)
	},
}

func init() {
	reformatCmd.AddCommand(fastaCmd)
}
