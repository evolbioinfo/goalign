package cmd

import (
	"errors"

	"github.com/fredericlemoine/goalign/io"
	"github.com/spf13/cobra"
)

var renameMap string
var renamerevert bool
var renameOutput string
var renameRegexp string
var renameReplace string

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
		var setregex, setreplace bool
		setregex = cmd.Flags().Changed("regexp")
		setreplace = cmd.Flags().Changed("replace")

		if renameMap == "none" {
			io.ExitWithMessage(errors.New("map file is not given"))
		}
		if setregex && !setreplace {
			io.ExitWithMessage(errors.New("--replace must be given with --regexp"))
		}

		// Read map file
		var namemap map[string]string
		var err error
		if !setregex {
			if namemap, err = readMapFile(renameMap, renamerevert); err != nil {
				io.ExitWithMessage(err)
			}
		} else {
			namemap = make(map[string]string)
		}

		aligns := readalign(infile)
		f := openWriteFile(renameOutput)
		for al := range aligns.Achan {
			if !setregex {
				al.Rename(namemap)
			} else {
				if err = al.RenameRegexp(renameRegexp, renameReplace, namemap); err != nil {
					io.ExitWithMessage(err)
				}
			}
			writeAlign(al, f)
		}
		f.Close()

		if setregex {
			writeNameMap(namemap, renameMap)
		}
	},
}

func init() {
	RootCmd.AddCommand(renameCmd)

	renameCmd.PersistentFlags().StringVarP(&renameMap, "map-file", "m", "none", "Name Mapping infile")
	renameCmd.PersistentFlags().BoolVarP(&renamerevert, "revert", "r", false, "Reverse orientation of mapfile")
	renameCmd.PersistentFlags().StringVarP(&renameOutput, "output", "o", "stdout", "renamed alignment output file")
	renameCmd.PersistentFlags().StringVarP(&renameRegexp, "regexp", "e", "none", "rename alignment using given regexp")
	renameCmd.PersistentFlags().StringVarP(&renameReplace, "replace", "b", "none", "replaces regexp matching strings by this string")
}
