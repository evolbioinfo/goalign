package cmd

import (
	"fmt"

	"github.com/evolbioinfo/goalign/align"
	"github.com/evolbioinfo/goalign/io"
	"github.com/spf13/cobra"
)

var statGapsFromStart bool
var statGapsFromEnd bool
var statGapsUnique bool

// charCmd represents the char command
var statGapsCmd = &cobra.Command{
	Use:   "gaps",
	Short: "Print gap stats on each alignment sequence",
	Long: `Print gap stats on each alignment sequence.

	By default, it prints, for each alignment sequence the number of gaps.

	- If --from-start is specified, then counts only gaps at sequence starts;
	- If --from-end is specified, then counts only gaps at sequence ends;
	- If --unique is specified, then counts only gaps that are unique in their column
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

		var numgapsunique []int
		if statGapsUnique {
			numgapsunique = al.NumGapsUniquePerSequence()
		}
		for i, s := range al.Sequences() {
			if statGapsFromStart {
				fmt.Printf("%s\t%d\n", s.Name(), s.NumGapsFromStart())
			} else if statGapsFromEnd {
				fmt.Printf("%s\t%d\n", s.Name(), s.NumGapsFromEnd())
			} else if statGapsUnique {
				fmt.Printf("%s\t%d\n", s.Name(), numgapsunique[i])
			} else {
				fmt.Printf("%s\t%d\n", s.Name(), s.NumGaps())
			}
		}

		return
	},
}

func init() {
	statGapsCmd.PersistentFlags().BoolVar(&statGapsFromStart, "from-start", false, "Count gaps in each sequence from start of sequences (until a non gap character is encountered)")
	statGapsCmd.PersistentFlags().BoolVar(&statGapsFromEnd, "from-end", false, "Count gaps in each sequence from end of sequences (until a non gap character is encountered)")
	statGapsCmd.PersistentFlags().BoolVar(&statGapsUnique, "unique", false, "Count, in each sequence, the number of gaps that are unique in a site")

	statsCmd.AddCommand(statGapsCmd)
}
