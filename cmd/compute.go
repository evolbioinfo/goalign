package cmd

import (
	"github.com/spf13/cobra"
)

var computeCmd = &cobra.Command{
	Use:   "compute",
	Short: "Different computations (distances, etc.)",
	Long: `Different computations (distances, etc.)
`,
}

func init() {
	RootCmd.AddCommand(computeCmd)
}
