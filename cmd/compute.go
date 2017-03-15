package cmd

import (
	"github.com/spf13/cobra"
)

var computeCmd = &cobra.Command{
	Use:   "compute",
	Short: "Different computations (distances, entropy, etc.)",
	Long: `Different computations (distances, entropy, etc.)
`,
}

func init() {
	RootCmd.AddCommand(computeCmd)
}
