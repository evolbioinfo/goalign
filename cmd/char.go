package cmd

import (
	"github.com/spf13/cobra"

	"github.com/fredericlemoine/goalign/io"
)

var charstatpersites bool

// charCmd represents the char command
var charCmd = &cobra.Command{
	Use:   "char",
	Short: "Prints frequence of different characters (aa/nt) of the alignment",
	Long: `Prints frequence of different characters (aa/nt) of the alignment.
May take a Phylip of Fasta input alignment.

Example of usages:

goalign stats char -i align.phylip -p
goalign stats char -i align.fasta
`,
	Run: func(cmd *cobra.Command, args []string) {
		al, _ := <-rootaligns.Achan
		if rootaligns.Err != nil {
			io.ExitWithMessage(rootaligns.Err)
		}
		if charstatpersites {
			printSiteCharStats(al)
		} else {
			printCharStats(al)
		}
	},
}

func init() {
	statsCmd.AddCommand(charCmd)
	charCmd.PersistentFlags().BoolVar(&charstatpersites, "per-sites", false, "Prints char statistics per alignment site")
}
