package cmd

import (
	"github.com/spf13/cobra"
)

// statsCmd represents the stats command
var statsCmd = &cobra.Command{
	Use:   "stats",
	Short: "Prints different characteristics of the alignment",
	Long: `Prints different characteristics of the alignment.

1 - Length
2 - Number of sequences
`,
}

func init() {
	RootCmd.AddCommand(statsCmd)
}
