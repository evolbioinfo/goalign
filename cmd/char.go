package cmd

import (
	"github.com/evolbioinfo/goalign/align"
	"github.com/evolbioinfo/goalign/io"
	"github.com/spf13/cobra"
)

var charstatpersites bool
var charstatpersequences bool

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
	RunE: func(cmd *cobra.Command, args []string) (err error) {
		var aligns *align.AlignChannel

		if aligns, err = readalign(infile); err != nil {
			io.LogError(err)
			return
		}

		for al := range aligns.Achan {
			if aligns.Err != nil {
				err = aligns.Err
				io.LogError(err)
				return
			}
			if charstatpersites {
				err = printSiteCharStats(al)
			} else if charstatpersequences {
				err = printSequenceCharStats(al)
			} else {
				printCharStats(al)
			}
		}
		return
	},
}

func init() {
	statsCmd.AddCommand(charCmd)
	charCmd.PersistentFlags().BoolVar(&charstatpersites, "per-sites", false, "Prints char statistics per alignment site (priority over --per-sequences)")
	charCmd.PersistentFlags().BoolVar(&charstatpersequences, "per-sequences", false, "Prints char statistics per alignment sequences")
}
