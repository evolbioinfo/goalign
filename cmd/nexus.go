package cmd

import (
	"github.com/evolbioinfo/goalign/align"
	"github.com/evolbioinfo/goalign/io"
	"github.com/evolbioinfo/goalign/io/utils"
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
			writeAlignNexus(al, f)
		}

		if aligns.Err != nil {
			err = aligns.Err
			io.LogError(err)
		}
		return
	},
}

func init() {
	reformatCmd.AddCommand(nexusCmd)
}
