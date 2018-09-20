package cmd

import (
	"fmt"

	"github.com/spf13/cobra"
)

var divideOutput string
var divideoutputFasta bool

// divideCmd represents the divide command
var divideCmd = &cobra.Command{
	Use:   "divide",
	Short: "Divide an input alignment in several output files",
	Long: `Divide an input alignment in several output files

If the alignment is in fasta format : will create 1 file
Otherwise, will create one file per alignment in the input file

-o : is the prefix of output files
if -o div, it will create files div_0.ph...div_n.ph

Output files will be in Phylip Format or in fasta format depending on -f

Example:

gotree divide -i align.ph -p -o out

`,
	Run: func(cmd *cobra.Command, args []string) {
		aligns := readalign(infile)
		i := 0
		for al := range aligns.Achan {
			if divideoutputFasta {
				f := openWriteFile(fmt.Sprintf("%s_%03d.fa", divideOutput, i))
				writeAlignFasta(al, f)
				f.Close()
			} else {
				f := openWriteFile(fmt.Sprintf("%s_%03d.ph", divideOutput, i))
				writeAlignPhylip(al, f)
				f.Close()
			}
			i++
		}
	},
}

func init() {
	RootCmd.AddCommand(divideCmd)
	divideCmd.PersistentFlags().StringVarP(&divideOutput, "output", "o", "prefix", "Divided alignment output files prefix")
	divideCmd.PersistentFlags().BoolVarP(&divideoutputFasta, "out-fasta", "f", false, "Output files in fasta format")

}
