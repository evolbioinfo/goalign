package cmd

import (
	"fmt"

	"github.com/fredericlemoine/goalign/align"
	"github.com/fredericlemoine/goalign/io"
	"github.com/spf13/cobra"
)

// charCmd represents the char command
var alphabetCmd = &cobra.Command{
	Use:   "alphabet",
	Short: "Prints the alphabet detected for the input alignment",
	Long: `Prints  the alphabet detected for the input alignment.

`,
	RunE: func(cmd *cobra.Command, args []string) (err error) {
		var aligns align.AlignChannel
		var seqs align.SeqBag

		if unaligned {
			if seqs, err = readsequences(infile); err != nil {
				io.LogError(err)
				return
			}
			fmt.Println(seqs.AlphabetStr())
		} else {

			if aligns, err = readalign(infile); err != nil {
				io.LogError(err)
				return
			} else {

				al, _ := <-aligns.Achan
				if aligns.Err != nil {
					err = aligns.Err
					io.LogError(err)
					return
				}
				fmt.Println(al.AlphabetStr())
			}
		}
		return
	},
}

func init() {
	alphabetCmd.PersistentFlags().BoolVar(&unaligned, "unaligned", false, "Considers sequences as unaligned and format fasta (phylip, nexus,... options are ignored)")
	statsCmd.AddCommand(alphabetCmd)
}
