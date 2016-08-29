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

1 - Length
2 - Number of sequences
`,
	Run: func(cmd *cobra.Command, args []string) {
		fmt.Fprintf(os.Stdout, "length\t%d\n", rootalign.Length())
		fmt.Fprintf(os.Stdout, "nseqs\t%d\n", rootalign.NbSequences())
		printCharStats(rootalign)
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
