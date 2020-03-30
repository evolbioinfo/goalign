package cmd

import (
	"fmt"

	"github.com/evolbioinfo/goalign/align"
	"github.com/evolbioinfo/goalign/io"
	"github.com/spf13/cobra"
)

var statMutationsRef string
var statMutationsUnique bool

// charCmd represents the char command
var statMutationsCmd = &cobra.Command{
	Use:   "mutations",
	Short: "Print mutations stats on each alignment sequence",
	Long: `Print mutations stats on each alignment sequence.

	- If --unique is specified, then counts only mutations (characters) that are unique in their column
	for the given sequence.

`,
	RunE: func(cmd *cobra.Command, args []string) (err error) {
		var aligns *align.AlignChannel

		if aligns, err = readalign(infile); err != nil {
			io.LogError(err)
			return
		}

		al, _ := <-aligns.Achan
		if aligns.Err != nil {
			err = aligns.Err
			io.LogError(err)
			return
		}

		var nummutationsunique []int
		if statMutationsUnique {
			nummutationsunique = al.NumMutationsUniquePerSequence()
		}
		for i, s := range al.Sequences() {
			if statMutationsUnique {
				fmt.Printf("%s\t%d\n", s.Name(), nummutationsunique[i])
			}
		}

		return
	},
}

func init() {
	statMutationsCmd.PersistentFlags().StringVar(&statMutationsRef, "ref-sequence", "none", "Count gaps in each sequence from end of sequences (until a non gap character is encountered)")
	statMutationsCmd.PersistentFlags().BoolVar(&statMutationsUnique, "unique", false, "Count, in each sequence, the number of mutations/characters that are unique in a site")

	statsCmd.AddCommand(statMutationsCmd)
}
