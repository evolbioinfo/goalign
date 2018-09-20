package cmd

import (
	"os"

	"github.com/fredericlemoine/goalign/align"
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

		a, _ := <-aligns.Achan
		if aligns.Err != nil {
			err = aligns.Err
			io.LogError(err)
			return
		}
		if reformatCleanNames {
			a.CleanNames()
		}
		writeAlignFasta(a, f)

		return
	},
}

func init() {
	reformatCmd.AddCommand(fastaCmd)
}
