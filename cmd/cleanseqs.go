package cmd

import (
	"fmt"

	"github.com/fredericlemoine/goalign/io"
	"github.com/spf13/cobra"
)

// cleanseqsCmd represents the cleanseqs command
var cleanseqsCmd = &cobra.Command{
	Use:   "seqs",
	Short: "Removes sequences with gaps",
	Long: `Removes sequences constituted of gaps

Removes sequences constitued of >= cutoff gap sites.

Exception for a cutoff of 0: removes sequencs constitued of > 0 gap sites.

Examples:
- With a cutoff of 0.5: a sequence with 5 gaps over 10 sites will be removed;
- With a cutoff of 0.5: a sequence with 4 gaps over 10 sites will not be removed;
- With a cutoff of 0.0 a site sequence 1 gap over 10 sites will be removed.

If cutoff is <0 or >1, it will be considered as 0, which means that every sequence with at least 1 gap
will be removed.`,
	Run: func(cmd *cobra.Command, args []string) {
		f := openWriteFile(cleanOutput)
		i := 0
		for al := range rootaligns {
			before := al.NbSequences()
			al.RemoveGapSeqs(cleanCutoff)
			after := al.NbSequences()
			writeAlign(al, f)
			if !cleanQuiet {
				io.PrintMessage(fmt.Sprintf("Alignment (%d) #seqs before cleaning=%d", i, before))
				io.PrintMessage(fmt.Sprintf("Alignment (%d) #seqs after cleaning=%d", i, after))
				io.PrintMessage(fmt.Sprintf("Alignment (%d) removed sequences=%d", i, before-after))
			}
		}
		f.Close()
	},
}

func init() {
	cleanCmd.AddCommand(cleanseqsCmd)
}
