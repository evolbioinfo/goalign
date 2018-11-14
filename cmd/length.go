package cmd

import (
	"fmt"

	"github.com/fredericlemoine/goalign/align"
	"github.com/fredericlemoine/goalign/io"
	"github.com/spf13/cobra"
)

// lengthCmd represents the length command
var lengthCmd = &cobra.Command{
	Use:   "length",
	Short: "Prints the length of sequences in the alignment",
	Long: `Prints the length of sequences in the alignment. 
May take a Phylip of Fasta input alignment.

If --unaligned is given, then length of all individual sequences is printed.

If the input alignment contains several alignments, will take all of them



Example of usages:

goalign stats length -i align.phylip -p
goalign stats length -i align.fasta

`,
	RunE: func(cmd *cobra.Command, args []string) (err error) {
		if unaligned {
			var seqs align.SeqBag

			if seqs, err = readsequences(infile); err != nil {
				io.LogError(err)
				return
			}

			seqs.IterateChar(func(name string, sequence []rune) {
				fmt.Println(name, "\t", len(sequence))
			})
		} else {
			var aligns align.AlignChannel

			if aligns, err = readalign(infile); err != nil {
				io.LogError(err)
				return
			}

			for al := range aligns.Achan {
				fmt.Println(al.Length())
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
	statsCmd.PersistentFlags().BoolVar(&unaligned, "unaligned", false, "Considers sequences as unaligned and format fasta (phylip, nexus,... options are ignored)")
	statsCmd.AddCommand(lengthCmd)
}
