package cmd

import (
	"github.com/spf13/cobra"
)

var cleanOutput string
var cleanCutoff float64
var cleanQuiet bool
var cleanChar string
var cleanIgnoreCase bool
var cleanIgnoreGaps bool
var cleanIgnoreNs bool

// cleanCmd represents the clean command
var cleanCmd = &cobra.Command{
	Use:   "clean",
	Short: "Removes gap sites or sequences",
	Long: `Removes sites or sequences constituted of gaps

Removes sites or sequences constitued of >= cutoff gap sites.

Exception for a cutoff of 0: removes sites constitued of > 0 gap sites.

Examples:
- With a cutoff of 0.5: a site with 5 gaps over 10 sequences will be removed;
- With a cutoff of 0.5: a site with 4 gaps over 10 sequences will not be removed;
- With a cutoff of 0.0 a site with 1 gap over 10 sequences will be removed.

If cutoff is <0 or >1, it will be considered as 0, which means that every site/sequence with at least 1 gap
will be removed.
`,
}

func init() {
	RootCmd.AddCommand(cleanCmd)
	cleanCmd.PersistentFlags().StringVarP(&cleanOutput, "output", "o", "stdout", "Cleaned alignment output file")
	cleanCmd.PersistentFlags().Float64VarP(&cleanCutoff, "cutoff", "c", 0, "Cutoff for deletion : 0 remove sites/sequences with > 0 given char, 1 remove sites/sequences with 100% given character)")
	cleanCmd.PersistentFlags().StringVar(&cleanChar, "char", "GAP", "The character the cutoff is applied to. May be GAP, MAJ, or any other character")
	cleanCmd.PersistentFlags().BoolVar(&cleanIgnoreCase, "ignore-case", false, "Ignore case of given character (--char) if non special character (GAP/-)")
	cleanCmd.PersistentFlags().BoolVar(&cleanIgnoreGaps, "ignore-gaps", false, "Ignore gaps for the % computation")
	cleanCmd.PersistentFlags().BoolVar(&cleanIgnoreNs, "ignore-n", false, "Ignore Ns for the % computation")
	cleanCmd.PersistentFlags().BoolVarP(&cleanQuiet, "quiet", "q", false, "Do not print results on stderr")
}
