package cmd

import (
	"github.com/spf13/cobra"
)

// reformatCmd represents the reformat command
var reformatCmd = &cobra.Command{
	Use:   "reformat",
	Short: "Reformats input alignment into phylip of fasta format",
	Long: `Reformats input alignment into phylip of fasta format.

Allows to convert formats between phylip and fasta, but also reformat 
fasta and phylip, such as 60 characters per line, etc.
`,
}

func init() {
	RootCmd.AddCommand(reformatCmd)
}
