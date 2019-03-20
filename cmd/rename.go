package cmd

import (
	"errors"
	"os"

	"github.com/fredericlemoine/goalign/align"
	"github.com/fredericlemoine/goalign/io"
	"github.com/spf13/cobra"
)

var renameMap string
var renamerevert bool
var renameCleanNames bool
var renameOutput string
var renameRegexp string
var renameReplace string

// renameCmd represents the rename command
var renameCmd = &cobra.Command{
	Use:   "rename",
	Short: "Rename sequences of the input alignment, given a map file",
	Long: `Rename sequences of the input tree, given a map file.

Different options:

1) A map file is given with --map-file: 
    It must be tab separated with the following columns:
      - Current name of the sequences
      - Desired new name of the sequences
    if --revert is given then it is the other way
    If a sequence name does not appear in the map file, it will
    not be renamed. 
    If a name that does not exist appears in the map file, it will
     not do anything.

2) --regexp and --replace are given:
    Sequence names are changed according to given regexp
    And mapping between old and new names is written in 
    the file potentially given with --map-file

3) --clean-names is given:
   Special characters in sequence names are replaced with '-'.
   Special characters are: \s\t[]();.,:|
   And mapping between old and new names is written in 
   the file potentially given with --map-file

In any case, option --unalign option will rename unaligned fasta files
while ignoring formatting options (phylip, etc.).


`,
	RunE: func(cmd *cobra.Command, args []string) (err error) {
		var setregex, setreplace bool
		var f *os.File
		var namemap map[string]string

		setregex = cmd.Flags().Changed("regexp")
		setreplace = cmd.Flags().Changed("replace")

		if setregex && !setreplace {
			err = errors.New("--replace must be given with --regexp")
			return
		}

		if f, err = openWriteFile(renameOutput); err != nil {
			io.LogError(err)
			return
		}
		defer closeWriteFile(f, renameOutput)

		// Read Map File
		if !setregex && !renameCleanNames {
			if renameMap == "none" {
				err = errors.New("map file is not given")
				return
			}
			if namemap, err = readMapFile(renameMap, renamerevert); err != nil {
				io.LogError(err)
				return
			}
		} else {
			namemap = make(map[string]string)
		}

		if unaligned {
			var seqs align.SeqBag

			if seqs, err = readsequences(infile); err != nil {
				io.LogError(err)
				return
			}
			if renameCleanNames {
				seqs.CleanNames(namemap)
			} else if setregex {
				if err = seqs.RenameRegexp(renameRegexp, renameReplace, namemap); err != nil {
					io.LogError(err)
					return
				}
			} else {
				seqs.Rename(namemap)
			}
			writeSequences(seqs, f)
		} else {
			var aligns *align.AlignChannel
			if aligns, err = readalign(infile); err != nil {
				io.LogError(err)
				return
			}
			for al := range aligns.Achan {
				if renameCleanNames {
					al.CleanNames(namemap)
				} else if setregex {
					if err = al.RenameRegexp(renameRegexp, renameReplace, namemap); err != nil {
						io.LogError(err)
						return
					}
				} else {
					al.Rename(namemap)
				}
				writeAlign(al, f)
			}
			if aligns.Err != nil {
				err = aligns.Err
				io.LogError(err)
			}
		}

		if (setregex || renameCleanNames) && renameMap != "None" {
			writeNameMap(namemap, renameMap)
		}

		return
	},
}

func init() {
	RootCmd.AddCommand(renameCmd)

	renameCmd.PersistentFlags().StringVarP(&renameMap, "map-file", "m", "none", "Name Mapping infile")
	renameCmd.PersistentFlags().BoolVarP(&renamerevert, "revert", "r", false, "Reverse orientation of mapfile")
	renameCmd.PersistentFlags().StringVarP(&renameOutput, "output", "o", "stdout", "renamed alignment output file")
	renameCmd.PersistentFlags().StringVarP(&renameRegexp, "regexp", "e", "none", "rename alignment using given regexp")
	renameCmd.PersistentFlags().StringVarP(&renameReplace, "replace", "b", "none", "replaces regexp matching strings by this string")
	renameCmd.PersistentFlags().BoolVar(&renameCleanNames, "clean-names", false, "Replaces special characters (tabs, spaces, newick characters) with '-' from input sequence names before writing output alignment")
	renameCmd.PersistentFlags().BoolVar(&unaligned, "unaligned", false, "Considers sequences as unaligned and format fasta (phylip, nexus,... options are ignored)")
}
