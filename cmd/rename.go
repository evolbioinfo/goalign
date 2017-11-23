package cmd

import (
	"errors"

	"github.com/fredericlemoine/goalign/io"
	"github.com/spf13/cobra"
)

var renameMap string
var renamerevert bool
var renameOutput string

// renameCmd represents the rename command
var renameCmd = &cobra.Command{
	Use:   "rename",
	Short: "Rename sequences of the input alignment, given a map file",
	Long: `Rename sequences of the input tree, given a map file.

Map file must be tab separated with columns:
1) Current name of the sequences
2) Desired new name of the sequences
(if --revert then it is the other way)

If a sequence name does not appear in the map file, it will not be renamed. 
If a name that does not exist appears in the map file, it will not do anything.
`,
	Run: func(cmd *cobra.Command, args []string) {

		if renameMap == "none" {
			io.ExitWithMessage(errors.New("map file is not given"))
		}

		// Read map file
		namemap, err := readMapFile(renameMap, renamerevert)
		if err != nil {
			io.ExitWithMessage(err)
		}

		f := openWriteFile(renameOutput)
		for al := range rootaligns.Achan {
			al.Rename(namemap)
			writeAlign(al, f)
		}
		f.Close()
	},
}

func init() {
	RootCmd.AddCommand(renameCmd)

	renameCmd.PersistentFlags().StringVarP(&renameMap, "map-file", "m", "none", "Name Mapping infile")
	renameCmd.PersistentFlags().BoolVarP(&renamerevert, "revert", "r", false, "Reverse orientation of mapfile")
	renameCmd.PersistentFlags().StringVarP(&renameOutput, "output", "o", "stdout", "renamed alignment output file")
}
