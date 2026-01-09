package cmd

import (
	"github.com/evolbioinfo/goalign/align"
	"github.com/evolbioinfo/goalign/io"
	"github.com/evolbioinfo/goalign/io/utils"
	"github.com/spf13/cobra"
)

// seqsCmd represents the seqs command
var seqsCmd = &cobra.Command{
	Use:   "seqs",
	Short: "Shuffles sequence order in alignment",
	Long: `Shuffle sequence order in alignment.

It may take a Fasta or Phylip alignment as input.

If the input alignment contains several alignments, will process all of them

Output a randomly reordered alignment. It does not
change the biological meaning of the alignment.

Example of usage:

goalign shuffle seqs -i align.phylip -p 
goalign shuffle seqs -i align.fasta

`,
	RunE: func(cmd *cobra.Command, args []string) (err error) {
		var f utils.StringWriterCloser

		if f, err = utils.OpenWriteFile(shuffleOutput); err != nil {
			io.LogError(err)
			return
		}
		defer utils.CloseWriteFile(f, shuffleOutput)

		if unaligned {
			var seqs align.SeqBag
			if seqs, err = readsequences(infile); err != nil {
				io.LogError(err)
				return
			}
			seqs.ShuffleSequences(globalRand)
			writeSequences(seqs, f)
		} else {
			var aligns *align.AlignChannel

			if aligns, err = readalign(infile); err != nil {
				io.LogError(err)
				return
			}

			for al := range aligns.Achan {
				al.ShuffleSequences(globalRand)
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
	seqsCmd.PersistentFlags().BoolVar(&unaligned, "unaligned", false, "Considers sequences as unaligned and format fasta (phylip, nexus,... options are ignored)")
	shuffleCmd.AddCommand(seqsCmd)
}
