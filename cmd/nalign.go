package cmd

import (
	"fmt"

	"github.com/spf13/cobra"
)

// nalignCmd represents the nalign command
var nalignCmd = &cobra.Command{
	Use:   "nalign",
	Short: "Prints the number of alignments in the input file",
	Long: `Prints the number of alignments in the input file

If the input file is in Fasta format, it should be 1
Otherwize, it may be > 1

Example:

goalign stats nalign -i align.ph -p

`,
	Run: func(cmd *cobra.Command, args []string) {
		naligns := 0
		for _ = range rootaligns {
			naligns++
		}
		fmt.Println(naligns)
	},
}

func init() {
	statsCmd.AddCommand(nalignCmd)
}
