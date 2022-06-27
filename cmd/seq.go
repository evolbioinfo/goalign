package cmd

import (
	"github.com/evolbioinfo/goalign/align"
	"github.com/evolbioinfo/goalign/io"
	"github.com/evolbioinfo/goalign/io/utils"
	"github.com/spf13/cobra"
)

var trimFromStart bool

// seqCmd represents the seq command
var seqCmd = &cobra.Command{
	Use:   "seq",
	Short: "Trims sequences of the alignment",
	Long: `Trims sequences of the alignemnt

It trims n (--nb-char, -n) characters from the beginning (--from-start, -s) or from the end (default) of the input alignment.

Example:
goalign trim seq -i align.fa -o trimed.fa -s -n 10

`,
	RunE: func(cmd *cobra.Command, args []string) (err error) {
		var aligns *align.AlignChannel
		var f utils.StringWriterCloser

		if aligns, err = readalign(infile); err != nil {
			io.LogError(err)
			return
		}
		if f, err = utils.OpenWriteFile(trimAlignOut); err != nil {
			io.LogError(err)
			return
		}
		defer utils.CloseWriteFile(f, trimAlignOut)

		for al := range aligns.Achan {
			if err = al.TrimSequences(trimNb, trimFromStart); err != nil {
				io.LogError(err)
				return
			} else {
				writeAlign(al, f)
			}
		}

		if aligns.Err != nil {
			err = aligns.Err
			io.LogError(err)
		}
		return
	},
}

func init() {
	trimCmd.AddCommand(seqCmd)
	seqCmd.PersistentFlags().IntVarP(&trimNb, "nb-char", "n", 1, "Number of characters to trim from sequences")
	seqCmd.PersistentFlags().BoolVarP(&trimFromStart, "from-start", "s", false, "If true: trims n char from the start, else from the end")
}
