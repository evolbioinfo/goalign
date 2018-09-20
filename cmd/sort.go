package cmd

import (
	"github.com/spf13/cobra"
)

var sortOutput string

// reformatCmd represents the reformat command
var sortCmd = &cobra.Command{
	Use:   "sort",
	Short: "sorts input alignment by sequence name",
	Long: `sorts input algignment by sequence name.
`,
	Run: func(cmd *cobra.Command, args []string) {
		aligns := readalign(infile)
		f := openWriteFile(sortOutput)
		for al := range aligns.Achan {
			al.Sort()
			writeAlign(al, f)
		}
		f.Close()
	},
}

func init() {
	RootCmd.AddCommand(sortCmd)
	sortCmd.PersistentFlags().StringVarP(&sortOutput, "output", "o", "stdout", "Sorted alignment output file")
}
