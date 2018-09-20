package cmd

import (
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
	Run: func(cmd *cobra.Command, args []string) {
		aligns := readalign(infile)
		f := openWriteFile(shuffleOutput)
		for al := range aligns.Achan {
			al.ShuffleSequences()
			writeAlign(al, f)
		}
		f.Close()
	},
}

func init() {
	shuffleCmd.AddCommand(seqsCmd)
}
