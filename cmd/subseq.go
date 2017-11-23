package cmd

import (
	"github.com/fredericlemoine/goalign/io"
	"github.com/spf13/cobra"
)

var subseqout string = "stdout"
var subseqstart int
var subseqlength int

// subseqCmd represents the subseq command
var subseqCmd = &cobra.Command{
	Use:   "subseq",
	Short: "Take a sub-alignment from the input alignment",
	Long: `Take a sub-alignment from the input alignment

It takes an alignment and extract sub-sequences from it, given
a start position (0-based inclusive) and a length.
If the length is after the end of the alignment, will stop at the 
end of the alignment.

For example:
goalign subseq -p -i al.phy -s 9 -e 10

This will extract a sub-alignment going from 10th position, with a length of 10.

The output format is the same than input format.
`,
	Run: func(cmd *cobra.Command, args []string) {
		out := openWriteFile(subseqout)
		for al := range rootaligns.Achan {
			subalign, err := al.SubAlign(subseqstart, subseqlength)
			if err != nil {
				io.ExitWithMessage(err)
			}
			writeAlign(subalign, out)
		}
		out.Close()
	},
}

func init() {
	RootCmd.AddCommand(subseqCmd)
	subseqCmd.PersistentFlags().StringVarP(&subseqout, "output", "o", "stdout", "Alignment output file")
	subseqCmd.PersistentFlags().IntVarP(&subseqstart, "start", "s", 0, "Start position (0-based inclusive)")
	subseqCmd.PersistentFlags().IntVarP(&subseqlength, "length", "l", 10, "Length of the sub alignment")
}
