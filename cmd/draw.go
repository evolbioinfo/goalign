package cmd

import (
	"github.com/spf13/cobra"
)

var drawOutput string

// drawCmd represents the draw command
var drawCmd = &cobra.Command{
	Use:   "draw",
	Short: "Draw alignments",
	Long:  `Draw alignments`,
}

func init() {
	RootCmd.AddCommand(drawCmd)
	drawCmd.PersistentFlags().StringVarP(&drawOutput, "output", "o", "stdout", "Alignment draw output file")
}
