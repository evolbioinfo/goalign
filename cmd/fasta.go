package cmd

import (
	"github.com/fredericlemoine/goalign/io"
	"github.com/spf13/cobra"
)

// fastaCmd represents the fasta command
var fastaCmd = &cobra.Command{
	Use:   "fasta",
	Short: "Reformats an input alignment into Fasta",
	Long: `Reformats an alignment into Fasta. 
It may take a Phylip of Fasta input alignment.

If the input alignment contains several alignments, will take the first one only


Example of usage:

goalign reformat fasta -i align.phylip -p
goalign reformat fasta -i align.fasta

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
		writeAlignFasta(a, f)
		f.Close()
	},
}

func init() {
	reformatCmd.AddCommand(fastaCmd)
}
