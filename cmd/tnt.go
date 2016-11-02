package cmd

import (
	"fmt"

	"github.com/spf13/cobra"
)

// tntCmd represents the tnt command
var tntCmd = &cobra.Command{
	Use:   "tnt",
	Short: "Reformats an input alignment into input data for TNT",
	Long: `Reformats an alignment into input data for TNT. 
It may take a Phylip of Fasta input alignment.

If the input alignment contains several alignments, will take the first one only

Example of usage:

goalign reformat tnt -i align.phylip -p
goalign reformat tnt -i align.fasta
`,
	Run: func(cmd *cobra.Command, args []string) {
		f := openWriteFile(reformatOutput)
		al := <-rootaligns
		f.WriteString("xread\n\n")
		f.WriteString("'Tnt input file'\n\n")
		f.WriteString(fmt.Sprintf("%d %d\n", al.Length(), al.NbSequences()))
		al.Iterate(func(name string, sequence string) {
			f.WriteString(fmt.Sprintf("%s %s\n", name, sequence))
		})
		f.WriteString(";\n")
		f.Close()
	},
}

func init() {
	reformatCmd.AddCommand(tntCmd)

}
