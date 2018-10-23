package cmd

import (
	"os"

	"github.com/fredericlemoine/goalign/align"
	"github.com/fredericlemoine/goalign/io"
	"github.com/spf13/cobra"
)

var translatePhase int
var translateOutput string

// translateCmd represents the addid command
var translateCmd = &cobra.Command{
	Use:   "translate",
	Short: "Translates an input alignment in amino acids",
	Long: `Translates an input alignment in amino acids.

If the input alignment is not nucleotides, then returns an error.

It is possible to drop a given number of characters from the start 
of the alignment, by specifying the '--phase' option.

It only translates using the standard genetic code so far.

`,
	RunE: func(cmd *cobra.Command, args []string) (err error) {
		var f *os.File

		if f, err = openWriteFile(translateOutput); err != nil {
			io.LogError(err)
			return
		}
		defer closeWriteFile(f, translateOutput)

		if unaligned {
			var seqs align.SeqBag

			if seqs, err = readsequences(infile); err != nil {
				io.LogError(err)
				return
			}
			if err = seqs.Translate(translatePhase); err != nil {
				io.LogError(err)
				return
			}
			writeSequences(seqs, f)
		} else {
			var aligns align.AlignChannel
			var al align.Alignment

			if aligns, err = readalign(infile); err != nil {
				io.LogError(err)
				return
			}
			for al = range aligns.Achan {
				if err = al.Translate(translatePhase); err != nil {
					io.LogError(err)
					return
				}
				writeAlign(al, f)
			}

			if aligns.Err != nil {
				err = aligns.Err
				io.LogError(err)
			}
		}

		return
	},
}

func init() {
	RootCmd.AddCommand(translateCmd)
	translateCmd.PersistentFlags().StringVarP(&translateOutput, "output", "o", "stdout", "Output translated alignment file")
	translateCmd.PersistentFlags().IntVar(&translatePhase, "phase", 0, "Number of characters to drop from the start of the alignment (if -1: Translate in the 3 phases, from positions 0, 1, and 2)")
	translateCmd.PersistentFlags().BoolVar(&unaligned, "unaligned", false, "Considers sequences as unaligned and format fasta (phylip, nexus,... options are ignored)")
}
