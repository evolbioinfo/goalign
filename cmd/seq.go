package cmd

import (
	"github.com/fredericlemoine/goalign/io"
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
	Run: func(cmd *cobra.Command, args []string) {
		aligns := readalign(infile)
		f := openWriteFile(trimAlignOut)
		for al := range aligns.Achan {
			if err := al.TrimSequences(trimNb, trimFromStart); err != nil {
				io.ExitWithMessage(err)
			} else {
				writeAlign(al, f)
			}
		}
		f.Close()
	},
}

func init() {
	trimCmd.AddCommand(seqCmd)
	seqCmd.PersistentFlags().IntVarP(&trimNb, "nb-char", "n", 1, "Number of characters to trim from sequences")
	seqCmd.PersistentFlags().BoolVarP(&trimFromStart, "from-start", "s", false, "If true: trims n char from the start, else from the end")
}
