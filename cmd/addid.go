package cmd

import (
	"os"

	"github.com/fredericlemoine/goalign/align"
	"github.com/fredericlemoine/goalign/io"
	"github.com/spf13/cobra"
)

var addIdOutput string
var addIdName string
var addIdRight bool

// addidCmd represents the addid command
var addidCmd = &cobra.Command{
	Use:   "addid",
	Short: "Adds a string to each sequence identifier of the input alignment",
	Long: `Adds a string to each sequence identifier of the input alignment.

The string may be added to the left or to the right of each sequence identifier.
`,
	RunE: func(cmd *cobra.Command, args []string) (err error) {
		var aligns align.AlignChannel
		var f *os.File

		if aligns, err = readalign(infile); err != nil {
			io.LogError(err)
			return
		}
		if f, err = openWriteFile(addIdOutput); err != nil {
			io.LogError(err)
			return
		}
		defer closeWriteFile(f, addIdOutput)

		for al := range aligns.Achan {
			al.AppendSeqIdentifier(addIdName, addIdRight)
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
	RootCmd.AddCommand(addidCmd)
	addidCmd.PersistentFlags().StringVarP(&addIdOutput, "out-align", "o", "stdout", "Renamed alignment output file")
	addidCmd.PersistentFlags().StringVarP(&addIdName, "name", "n", "none", "String to add to sequence names")
	addidCmd.PersistentFlags().BoolVarP(&addIdRight, "right", "r", false, "Adds the String on the right of sequence names (otherwise, adds to left)")
}
