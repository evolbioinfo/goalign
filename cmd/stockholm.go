package cmd

import (
	"github.com/evolbioinfo/goalign/align"
	"github.com/evolbioinfo/goalign/io"
	"github.com/evolbioinfo/goalign/io/utils"
	"github.com/spf13/cobra"
)

var stockholmCmd = &cobra.Command{
	Use:   "stockholm",
	Short: "Reformats an input alignment into Stockholm format",
	Long: `Reformats an alignment into Stockholm format. 
If the input alignment contains several alignments, will take the first one only


Example of usage:

goalign reformat stockholm -i align.phylip -p
goalign reformat stockholm -i align.fasta

`,
	RunE: func(cmd *cobra.Command, args []string) (err error) {
		var f utils.StringWriterCloser

		if f, err = utils.OpenWriteFile(reformatOutput); err != nil {
			io.LogError(err)
			return
		}
		defer utils.CloseWriteFile(f, reformatOutput)

		var aligns *align.AlignChannel

		if aligns, err = readalign(infile); err != nil {
			io.LogError(err)
			return
		}

		a := <-aligns.Achan
		if aligns.Err != nil {
			err = aligns.Err
			io.LogError(err)
			return
		}
		if reformatCleanNames {
			a.CleanNames(nil)
		}
		writeAlignStockholm(a, f)
		return
	},
}

func init() {
	reformatCmd.AddCommand(stockholmCmd)
}
