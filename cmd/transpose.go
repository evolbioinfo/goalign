package cmd

import (
	"github.com/evolbioinfo/goalign/align"
	"github.com/evolbioinfo/goalign/io"
	"github.com/evolbioinfo/goalign/io/utils"
	"github.com/spf13/cobra"
)

var transposeOutput string

// transposeCmd represents the addid command
var transposeCmd = &cobra.Command{
	Use:   "transpose",
	Short: "Transpose an input alignment",
	Long: `Transposes an input alignment such that the sequences 
	become the sites and the sites become the sequence.
`,
	RunE: func(cmd *cobra.Command, args []string) (err error) {
		var f utils.StringWriterCloser
		var aligns *align.AlignChannel
		var al, tr align.Alignment

		if f, err = utils.OpenWriteFile(transposeOutput); err != nil {
			io.LogError(err)
			return
		}
		defer utils.CloseWriteFile(f, transposeOutput)

		if aligns, err = readalign(infile); err != nil {
			io.LogError(err)
			return
		}
		for al = range aligns.Achan {
			if tr, err = al.Transpose(); err != nil {
				io.LogError(err)
				return
			}
			writeAlign(tr, f)
		}

		if aligns.Err != nil {
			err = aligns.Err
			io.LogError(err)
		}

		return
	},
}

func init() {
	RootCmd.AddCommand(transposeCmd)
	transposeCmd.PersistentFlags().StringVarP(&transposeOutput, "output", "o", "stdout", "Output transposed alignment")
}
