package cmd

import (
	"github.com/fredericlemoine/goalign/io"
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
	Run: func(cmd *cobra.Command, args []string) {
		aligns := readalign(infile)
		f := openWriteFile(reformatOutput)
		a, _ := <-aligns.Achan
		if aligns.Err != nil {
			io.ExitWithMessage(aligns.Err)
		}
		if reformatCleanNames {
			a.CleanNames()
		}
		writeAlignClustal(a, f)
		f.Close()
	},
}

func init() {
	reformatCmd.AddCommand(clustalCmd)
}
