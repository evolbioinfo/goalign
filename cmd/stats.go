package cmd

import (
	"fmt"
	"github.com/fredericlemoine/goalign/align"
	"github.com/spf13/cobra"
	"os"
	"sort"
)

// statsCmd represents the stats command
var statsCmd = &cobra.Command{
	Use:   "stats",
	Short: "Prints different characteristics of the alignment",
	Long: `Prints different characteristics of the alignment.

1. Length of alignment;
2. Number of sequences;
3. Average number of alleles per site;
4. Number of variables sites (does ot take into account gaps or special characters);
5. Character frequencies.

If the input alignment contains several alignments, will process all of them

`,
	Run: func(cmd *cobra.Command, args []string) {
		for al := range rootaligns.Achan {
			fmt.Fprintf(os.Stdout, "length\t%d\n", al.Length())
			fmt.Fprintf(os.Stdout, "nseqs\t%d\n", al.NbSequences())
			fmt.Fprintf(os.Stdout, "avgalleles\t%.4f\n", al.AvgAllelesPerSite())
			fmt.Fprintf(os.Stdout, "variable sites\t%d\n", al.NbVariableSites())
			printCharStats(al)
		}
	},
}

func printCharStats(rootalign align.Alignment) {
	charmap := rootalign.CharStats()
	keys := make([]string, 0, len(charmap))
	var total int64 = 0
	for k, v := range charmap {
		keys = append(keys, string(k))
		total += v
	}
	sort.Strings(keys)

	fmt.Fprintf(os.Stdout, "char\tnb\tfreq\n")
	for _, k := range keys {
		nb := charmap[rune(k[0])]
		fmt.Fprintf(os.Stdout, "%s\t%d\t%f\n", k, nb, float64(nb)/float64(total))
	}
}

func init() {
	RootCmd.AddCommand(statsCmd)
}
