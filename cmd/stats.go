package cmd

import (
	"fmt"
	"os"
	"sort"

	"github.com/fredericlemoine/goalign/align"
	"github.com/fredericlemoine/goalign/io"
	"github.com/spf13/cobra"
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
	RunE: func(cmd *cobra.Command, args []string) (err error) {
		var aligns align.AlignChannel

		if aligns, err = readalign(infile); err != nil {
			io.LogError(err)
			return
		}

		for al := range aligns.Achan {
			fmt.Fprintf(os.Stdout, "length\t%d\n", al.Length())
			fmt.Fprintf(os.Stdout, "nseqs\t%d\n", al.NbSequences())
			fmt.Fprintf(os.Stdout, "avgalleles\t%.4f\n", al.AvgAllelesPerSite())
			fmt.Fprintf(os.Stdout, "variable sites\t%d\n", al.NbVariableSites())
			printCharStats(al)
		}

		if aligns.Err != nil {
			err = aligns.Err
			io.LogError(err)
		}
		return
	},
}

func printCharStats(align align.Alignment) {
	charmap := align.CharStats()
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

func printSiteCharStats(align align.Alignment) (err error) {
	var sitemap map[rune]int

	charmap := align.CharStats()
	keys := make([]string, 0, len(charmap))
	for k, _ := range charmap {
		keys = append(keys, string(k))
	}
	sort.Strings(keys)
	fmt.Fprintf(os.Stdout, "site")
	for _, v := range keys {
		fmt.Fprintf(os.Stdout, "\t%s", v)
	}
	fmt.Fprintf(os.Stdout, "\n")
	for site := 0; site < align.Length(); site++ {
		if sitemap, err = align.CharStatsSite(site); err != nil {
			return
		}
		fmt.Fprintf(os.Stdout, "%d", site)
		for _, k := range keys {
			nb := sitemap[rune(k[0])]
			fmt.Fprintf(os.Stdout, "\t%d", nb)
		}
		fmt.Fprintf(os.Stdout, "\n")
	}
	return
}

// Prints the Character with the most frequency
// for each site of the alignment
func printMaxCharStats(align align.Alignment) {
	maxchars, occur := align.MaxCharStats()

	fmt.Fprintf(os.Stdout, "site\tchar\tnb\n")
	for i, c := range maxchars {
		fmt.Fprintf(os.Stdout, "%d\t%c\t%d\n", i, c, occur[i])
	}
}

func init() {
	RootCmd.AddCommand(statsCmd)
}
