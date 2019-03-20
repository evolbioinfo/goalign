package cmd

import (
	"fmt"
	"os"

	"github.com/fredericlemoine/goalign/align"
	"github.com/fredericlemoine/goalign/io"
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
	RunE: func(cmd *cobra.Command, args []string) (err error) {
		var aligns *align.AlignChannel
		var f *os.File

		if aligns, err = readalign(infile); err != nil {
			io.LogError(err)
			return
		}

		i := 0
		for al := range aligns.Achan {
			if divideoutputFasta {
				if f, err = openWriteFile(fmt.Sprintf("%s_%03d.fa", divideOutput, i)); err != nil {
					io.LogError(err)
					return
				}
				writeAlignFasta(al, f)
				f.Close()
			} else {
				if f, err = openWriteFile(fmt.Sprintf("%s_%03d.ph", divideOutput, i)); err != nil {
					io.LogError(err)
					return
				}
				writeAlignPhylip(al, f)
				f.Close()
			}
			i++
		}

		if aligns.Err != nil {
			err = aligns.Err
			io.LogError(err)
		}
		return
	},
}

func init() {
	RootCmd.AddCommand(divideCmd)
	divideCmd.PersistentFlags().StringVarP(&divideOutput, "output", "o", "prefix", "Divided alignment output files prefix")
	divideCmd.PersistentFlags().BoolVarP(&divideoutputFasta, "out-fasta", "f", false, "Output files in fasta format")

}
