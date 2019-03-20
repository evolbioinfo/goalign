package cmd

import (
	"fmt"

	"github.com/fredericlemoine/goalign/align"
	"github.com/fredericlemoine/goalign/io"
	"github.com/spf13/cobra"
)

// taxaCmd represents the taxa command
var taxaCmd = &cobra.Command{
	Use:   "taxa",
	Short: "Prints index (position) and name of taxa of the alignment file",
	Long: `Prints index (position) and name of taxa of the alignment file.
May take a Phylip of Fasta input alignment.

If the input alignment contains several alignments, will take only the first one

Example of usages:

goalign stats taxa -i align.phylip -p
goalign stats taxa -i align.fasta
`,
	RunE: func(cmd *cobra.Command, args []string) (err error) {

		if unaligned {
			var seqs align.SeqBag

			if seqs, err = readsequences(infile); err != nil {
				io.LogError(err)
				return
			}
			i := 0
			seqs.Iterate(func(name string, sequence string) {
				fmt.Print(fmt.Sprintf("%d\t%s\n", i, name))
				i++
			})
		} else {
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

			i := 0
			al.Iterate(func(name string, sequence string) {
				fmt.Print(fmt.Sprintf("%d\t%s\n", i, name))
				i++
			})
		}
		return
	},
}

func init() {
	statsCmd.PersistentFlags().BoolVar(&unaligned, "unaligned", false, "Considers sequences as unaligned and format fasta (phylip, nexus,... options are ignored)")
	statsCmd.AddCommand(taxaCmd)
}
