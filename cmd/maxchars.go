package cmd

import (
	"github.com/evolbioinfo/goalign/align"
	"github.com/evolbioinfo/goalign/io"
	"github.com/spf13/cobra"
)

var maxCharExcludeGaps bool
var maxCharIgnoreGaps bool
var maxCharIgnoreNs bool

// charCmd represents the char command
var maxCharCmd = &cobra.Command{
	Use:   "maxchar",
	Short: "Prints the character with the highest occcurence for each site of the alignment",
	Long: `Prints the character with the highest occcurence for each site of the alignment.

Ouput format: Tabulated with columns:
1) Site index (0...)
2) Character with maximum occurence
3) Number of occurence of this character

Example of usages:

goalign stats maxchar -i align.phylip -p
goalign stats maxchar -i align.fasta
`,
	RunE: func(cmd *cobra.Command, args []string) (err error) {
		var aligns *align.AlignChannel
		maxCharIgnoreGaps = maxCharIgnoreGaps || maxCharExcludeGaps

		if aligns, err = readalign(infile); err != nil {
			io.LogError(err)
			return
		}

		al := <-aligns.Achan
		if aligns.Err != nil {
			err = aligns.Err
			io.LogError(err)
			return
		}
		printMaxCharStats(al, maxCharIgnoreGaps, maxCharIgnoreNs)

		return
	},
}

func init() {
	statsCmd.AddCommand(maxCharCmd)

	maxCharCmd.PersistentFlags().BoolVar(&maxCharExcludeGaps, "exclude-gaps", false, "Ignore gaps in the majority computation (for backward compatibility, will be removed in future releases)")
	maxCharCmd.PersistentFlags().BoolVar(&maxCharIgnoreGaps, "ignore-gaps", false, "Ignore gaps in the majority computation")
	maxCharCmd.PersistentFlags().BoolVar(&maxCharIgnoreNs, "ignore-n", false, "Ignore Ns in the majority computation")
}
