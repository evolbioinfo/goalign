package cmd

import (
	"github.com/evolbioinfo/goalign/align"
	"github.com/evolbioinfo/goalign/io"
	"github.com/evolbioinfo/goalign/io/utils"
	"github.com/spf13/cobra"
)

// clustalCmd : to reformat in clustal format
var clustalCmd = &cobra.Command{
	Use:   "clustal",
	Short: "Reformats an input alignment into Clustal format",
	Long: `Reformats an alignment into Clustal format. 
It may take a Phylip, Fasta, Nexus, or Clustal input alignment.

Example of usage:

goalign reformat clustal -i align.phylip -p
goalign reformat clustal -i align.fasta

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

		a := <-aligns.Achan
		if aligns.Err != nil {
			err = aligns.Err
			io.LogError(err)
			return
		}
		if reformatCleanNames {
			a.CleanNames(nil)
		}
		writeAlignClustal(a, f)

		return
	},
}

func init() {
	reformatCmd.AddCommand(clustalCmd)
}
