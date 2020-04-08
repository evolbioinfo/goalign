package cmd

import (
	"fmt"
	"os"

	"github.com/evolbioinfo/goalign/align"
	"github.com/evolbioinfo/goalign/io"
	"github.com/spf13/cobra"
)

// tntCmd represents the tnt command
var tntCmd = &cobra.Command{
	Use:   "tnt",
	Short: "Reformats an input alignment into input data for TNT",
	Long: `Reformats an alignment into input data for TNT. 
It may take a Phylip of Fasta input alignment.

If the input alignment contains several alignments, will take the first one only

Example of usage:

goalign reformat tnt -i align.phylip -p
goalign reformat tnt -i align.fasta
`,
	RunE: func(cmd *cobra.Command, args []string) (err error) {
		var aligns *align.AlignChannel
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

		al, _ := <-aligns.Achan
		if aligns.Err != nil {
			err = aligns.Err
			io.LogError(err)
			return
		}

		if reformatCleanNames {
			al.CleanNames(nil)
		}
		f.WriteString("xread\n\n")
		f.WriteString("'Tnt input file'\n\n")
		f.WriteString(fmt.Sprintf("%d %d\n", al.Length(), al.NbSequences()))
		al.Iterate(func(name string, sequence string) bool {
			f.WriteString(fmt.Sprintf("%s %s\n", name, sequence))
			return false
		})
		f.WriteString(";\n")

		return
	},
}

func init() {
	reformatCmd.AddCommand(tntCmd)

}
