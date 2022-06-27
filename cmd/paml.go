package cmd

import (
	"github.com/evolbioinfo/goalign/align"
	"github.com/evolbioinfo/goalign/io"
	"github.com/evolbioinfo/goalign/io/utils"
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
		var aligns *align.AlignChannel
		var f utils.StringWriterCloser

		if aligns, err = readalign(infile); err != nil {
			io.LogError(err)
			return
		}
		if f, err = utils.OpenWriteFile(reformatOutput); err != nil {
			io.LogError(err)
			return
		}
		defer utils.CloseWriteFile(f, reformatOutput)

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
