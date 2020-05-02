package cmd

import (
	"fmt"
	"os"

	"github.com/evolbioinfo/goalign/align"
	"github.com/evolbioinfo/goalign/io"
	"github.com/spf13/cobra"
)

var cleanEnds bool
var cleanChar string

// cleansitesCmd represents the cleansites command
var cleansitesCmd = &cobra.Command{
	Use:   "sites",
	Short: "Removes sites with specific characters",
	Long: `Removes sites constituted of specific characters

Removes sites constitued of >= cutoff specific characters. This characters can be :

1. Gap (--char=GAP or --char=-, default)
2. Any other character X specified by --char=X (case sensitive)
3. The most abundant character in the site --char=MAJ (including gaps)

Exception for a cutoff of 0: removes sites constitued of > 0 specified character (with --char=MAJ, then will remove all columns).

Examples:
- With a cutoff of 0.5: a site with 5 specified characters over 10 sequences will be removed;
- With a cutoff of 0.5: a site with 4 specified characters over 10 sequences will not be removed;
- With a cutoff of 0.0 a site with 1 specified over 10 sequences will be removed.

If cutoff is <0 or >1, it will be considered as 0, which means that every site with at least 1 specified character
will be removed.`,
	RunE: func(cmd *cobra.Command, args []string) (err error) {
		var aligns *align.AlignChannel
		var nbstart, nbend int
		var f *os.File

		if aligns, err = readalign(infile); err != nil {
			io.LogError(err)
			return
		}
		if f, err = openWriteFile(cleanOutput); err != nil {
			io.LogError(err)
			return
		}
		defer closeWriteFile(f, cleanOutput)

		i := 0
		char := ""

		for al := range aligns.Achan {
			beforelength := al.Length()

			if cleanChar == string(align.GAP) || cleanChar == "GAP" {
				char = "gaps"
				nbstart, nbend = al.RemoveGapSites(cleanCutoff, cleanEnds)
			} else if cleanChar == "MAJ" {
				char = "maj"
				nbstart, nbend = al.RemoveMajorityCharacterSites(cleanCutoff, cleanEnds)
			} else {
				//single character
				c := []rune(cleanChar)
				if len(c) != 1 {
					err = fmt.Errorf("--char should be a single character")
					io.LogError(err)
					return
				}
				char = string(c[0])
				nbstart, nbend = al.RemoveCharacterSites(c[0], cleanCutoff, cleanEnds)
			}
			afterlength := al.Length()
			writeAlign(al, f)
			if !cleanQuiet {
				io.PrintSimpleMessage(fmt.Sprintf("Alignment (%d) length before cleaning=%d", i, beforelength))
				io.PrintSimpleMessage(fmt.Sprintf("Alignment (%d) length after cleaning=%d", i, afterlength))
				io.PrintSimpleMessage(fmt.Sprintf("Alignment (%d) number of %s=%d", i, char, beforelength-afterlength))
				io.PrintSimpleMessage(fmt.Sprintf("Alignment (%d) number of start %s=%d", i, char, nbstart))
				io.PrintSimpleMessage(fmt.Sprintf("Alignment (%d) number of end %s=%d", i, char, nbend))
			}
		}

		if aligns.Err != nil {
			err = aligns.Err
			io.LogError(err)
		}
		return
	},
}

func init() {
	cleansitesCmd.PersistentFlags().BoolVar(&cleanEnds, "ends", false, "If true, then only remove consecutive gap positions from alignment start and end")
	cleansitesCmd.PersistentFlags().StringVar(&cleanChar, "char", "GAP", "The character the cutoff is applied to. May be GAP, MAJ, or any other character")

	cleanCmd.AddCommand(cleansitesCmd)
}
