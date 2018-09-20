package cmd

import (
	"os"

	"github.com/fredericlemoine/goalign/align"
	"github.com/fredericlemoine/goalign/io"
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
				al.CleanNames()
			}
			writeAlignPhylip(al, f)
		}

		if aligns.Err != nil {
			err = aligns.Err
			io.LogError(err)
		}
		return
	},
}

func init() {
	reformatCmd.AddCommand(phylipCmd)
}
