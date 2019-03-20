package cmd

import (
	"os"

	"github.com/fredericlemoine/goalign/align"
	"github.com/fredericlemoine/goalign/io"
	"github.com/spf13/cobra"
)

var sortOutput string

// reformatCmd represents the reformat command
var sortCmd = &cobra.Command{
	Use:   "sort",
	Short: "sorts input alignment by sequence name",
	Long: `sorts input algignment by sequence name.
`,
	RunE: func(cmd *cobra.Command, args []string) (err error) {
		var aligns *align.AlignChannel
		var f *os.File

		if aligns, err = readalign(infile); err != nil {
			io.LogError(err)
			return
		}
		if f, err = openWriteFile(sortOutput); err != nil {
			io.LogError(err)
			return
		}
		defer closeWriteFile(f, sortOutput)

		for al := range aligns.Achan {
			al.Sort()
			writeAlign(al, f)
		}

		if aligns.Err != nil {
			err = aligns.Err
			io.LogError(err)
		}
		return
	},
}

func init() {
	RootCmd.AddCommand(sortCmd)
	sortCmd.PersistentFlags().StringVarP(&sortOutput, "output", "o", "stdout", "Sorted alignment output file")
}
