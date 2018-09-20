package cmd

import (
	"github.com/spf13/cobra"
)

// tntCmd represents the tnt command
var pamlCmd = &cobra.Command{
	Use:   "paml",
	Short: "Reformats an input alignment into input data for PAML",
	Long: `Reformats an alignment into input data for PAML. 
It may take a Phylip of Fasta input alignment.

If the input alignment contains several alignments, will take the first one only

Example of usage:

goalign reformat paml -i align.phylip -p
goalign reformat paml -i align.fasta
`,
	Run: func(cmd *cobra.Command, args []string) {
		aligns := readalign(infile)
		f := openWriteFile(reformatOutput)
		for al := range aligns.Achan {
			if reformatCleanNames {
				al.CleanNames()
			}
			writeAlignPaml(al, f)
		}
		f.Close()
	},
}

func init() {
	reformatCmd.AddCommand(pamlCmd)

}
