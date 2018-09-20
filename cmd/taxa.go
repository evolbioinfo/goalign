package cmd

import (
	"fmt"

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
	Run: func(cmd *cobra.Command, args []string) {
		aligns := readalign(infile)
		al, _ := <-aligns.Achan
		if aligns.Err != nil {
			io.ExitWithMessage(aligns.Err)
		}

		i := 0
		al.Iterate(func(name string, sequence string) {
			fmt.Print(fmt.Sprintf("%d\t%s\n", i, name))
			i++
		})
	},
}

func init() {
	statsCmd.AddCommand(taxaCmd)
}
