package cmd

import (
	"github.com/spf13/cobra"
)

// nexusCmd represents the nexus command
var nexusCmd = &cobra.Command{
	Use:   "nexus",
	Short: "Reformats an input alignment into nexus",
	Long: `Reformats an alignment into nexus format. 
It may take a Phylip of Fasta input alignment.

If the input alignment contains several alignments, will take all of them

Example of usage:

goalign reformat nexus -i align.phylip -p
goalign reformat nexus -i align.fasta

`,
	Run: func(cmd *cobra.Command, args []string) {
		f := openWriteFile(reformatOutput)
		for al := range rootaligns.Achan {
			//fmt.Println("ALIGN" + fmt.Sprintf("%d", al.NbSequences()))
			writeAlignNexus(al, f)
		}
		f.Close()
	},
}

func init() {
	reformatCmd.AddCommand(nexusCmd)
}
