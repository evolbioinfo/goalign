package cmd

import (
	"os"

	"github.com/fredericlemoine/goalign/align"
	"github.com/fredericlemoine/goalign/io"
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
	RunE: func(cmd *cobra.Command, args []string) (err error) {
		var aligns align.AlignChannel
		var f *os.File

		if aligns, err = readalign(infile); err != nil {
			io.LogError(err)
			return
		}
		if f, err = openWriteFile(reformatOutput); err != nil {
			io.LogError(err)
			return
		}
		defer closeWriteFile(f, reformatOutput)

		for al := range aligns.Achan {
			if reformatCleanNames {
				al.CleanNames(nil)
			}
			writeAlignPaml(al, f)
		}

		if aligns.Err != nil {
			err = aligns.Err
			io.LogError(err)
		}
		return
	},
}

func init() {
	reformatCmd.AddCommand(pamlCmd)

}
