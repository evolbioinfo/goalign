package cmd

import (
	"fmt"

	"github.com/fredericlemoine/goalign/align"
	"github.com/fredericlemoine/goalign/io"
	"github.com/spf13/cobra"
)

// nseqCmd represents the nseq command
var nseqCmd = &cobra.Command{
	Use:   "nseq",
	Short: "Prints the number of sequences in the alignment",
	Long: `Prints the number of sequences in the alignment. 
May take a Phylip of Fasta input alignment.

If the input alignment contains several alignments, will process all of them

Example of usages:

goalign stats nseq -i align.phylip -p
goalign stats nseq -i align.fasta
`,
	RunE: func(cmd *cobra.Command, args []string) (err error) {
		var aligns *align.AlignChannel

		if aligns, err = readalign(infile); err != nil {
			io.LogError(err)
			return
		}

		for al := range aligns.Achan {
			fmt.Println(al.NbSequences())
		}

		if aligns.Err != nil {
			err = aligns.Err
			io.LogError(err)
		}
		return
	},
}

func init() {
	statsCmd.AddCommand(nseqCmd)
}
