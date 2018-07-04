package cmd

import (
	"fmt"

	"github.com/fredericlemoine/goalign/version"
	"github.com/spf13/cobra"
)

// versionCmd represents the version command
var versionCmd = &cobra.Command{
	Use:   "version",
	Short: "Prints the current version of goalign",
	Long:  `Prints the current version of goalign`,
	PersistentPreRun: func(cmd *cobra.Command, args []string) {
	},
	Run: func(cmd *cobra.Command, args []string) {
		fmt.Println(version.Version)
	},
}

func init() {
	RootCmd.AddCommand(versionCmd)
}
