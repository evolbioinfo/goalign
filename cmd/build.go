package cmd

import (
	"github.com/spf13/cobra"
)

// buildCmd represents the build command
var buildCmd = &cobra.Command{
	Use:   "build",
	Short: "Command to build output files : bootstrap for example",
	Long:  `Only one possibility so far: building bootstrap alignment.`,
}

func init() {
	RootCmd.AddCommand(buildCmd)
}
