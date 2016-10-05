package cmd

import (
	"github.com/spf13/cobra"
)

var phylipstrict bool

// phylipCmd represents the phylip command
var phylipCmd = &cobra.Command{
	Use:   "phylip",
	Short: "Reformats an input alignment into Phylip",
	Long: `Reformats an alignment into Phylip. 
It may take a Phylip of Fasta input alignment.

If the input alignment contains several alignments, will take all of them

Example of usage:

goalign reformat phylip -i align.phylip -p
goalign reformat phylip -i align.fasta

`,
	Run: func(cmd *cobra.Command, args []string) {
		f := openWriteFile(reformatOutput)
		for al := range rootaligns {
			//fmt.Println("ALIGN" + fmt.Sprintf("%d", al.NbSequences()))
			writeAlignPhylip(al, f, phylipstrict)
		}
		f.Close()
	},
}

func init() {
	reformatCmd.AddCommand(phylipCmd)
	reformatCmd.PersistentFlags().BoolVarP(&phylipstrict, "strict", "s", false, "If it is strict phylip format")
}
