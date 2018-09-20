package cmd

import (
	"fmt"
	"github.com/spf13/cobra"

	"github.com/fredericlemoine/goalign/io"
)

// cleansitesCmd represents the cleansites command
var cleansitesCmd = &cobra.Command{
	Use:   "sites",
	Short: "Removes sites with gaps",
	Long: `Removes sites constituted of gaps

Removes sites constitued of >= cutoff gap sites.

Exception for a cutoff of 0: removes sites constitued of > 0 gap sites.

Examples:
- With a cutoff of 0.5: a site with 5 gaps over 10 sequences will be removed;
- With a cutoff of 0.5: a site with 4 gaps over 10 sequences will not be removed;
- With a cutoff of 0.0 a site with 1 gap over 10 sequences will be removed.

If cutoff is <0 or >1, it will be considered as 0, which means that every site with at least 1 gap
will be removed.`,
	Run: func(cmd *cobra.Command, args []string) {
		aligns := readalign(infile)
		f := openWriteFile(cleanOutput)
		i := 0
		for al := range aligns.Achan {
			beforelength := al.Length()
			al.RemoveGapSites(cleanCutoff)
			afterlength := al.Length()
			writeAlign(al, f)
			if !cleanQuiet {
				io.PrintMessage(fmt.Sprintf("Alignment (%d) length before cleaning=%d", i, beforelength))
				io.PrintMessage(fmt.Sprintf("Alignment (%d) length after cleaning=%d", i, afterlength))
				io.PrintMessage(fmt.Sprintf("Alignment (%d) number of gaps=%d", i, beforelength-afterlength))
			}
		}
		f.Close()
	},
}

func init() {
	cleanCmd.AddCommand(cleansitesCmd)
}
