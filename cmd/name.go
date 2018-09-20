package cmd

import (
	"os"

	"github.com/fredericlemoine/goalign/align"
	"github.com/fredericlemoine/goalign/io"
	"github.com/spf13/cobra"
)

var trimMapout string
var trimAuto bool

// nameCmd represents the name command
var nameCmd = &cobra.Command{
	Use:   "name",
	Short: "This command trims names of sequences",
	Long: `This command trims names of sequences.

If the input alignment contains several alignments, will process only the first one

It trims sequence names to n characters. It will output mapping between 
old names and new names into a map file as well as the new alignment.
If n is > than seq name, it will add 000 at the end.
If n is < than seq name, it will trim. Names may be identical after this step.
In this case, it will add a unique identifier to the identical names.

Example of usage:

goalign trim name -i align.phy -p -n 10 -m map.txt

Id -a is given, then names are generated with the pattern "S000<i>".
`,
	RunE: func(cmd *cobra.Command, args []string) (err error) {
		var aligns align.AlignChannel
		var f *os.File

		if aligns, err = readalign(infile); err != nil {
			io.LogError(err)
			return
		}
		if f, err = openWriteFile(trimAlignOut); err != nil {
			io.LogError(err)
			return
		}
		defer closeWriteFile(f, trimAlignOut)

		namemap := make(map[string]string)
		curid := 1
		for al := range aligns.Achan {
			if aligns.Err != nil {
				err = aligns.Err
				io.LogError(err)
				return
			}

			if trimAuto {
				if err = al.TrimNamesAuto(namemap, &curid); err != nil {
					io.LogError(err)
					return
				}
			} else {
				if err = al.TrimNames(namemap, trimNb); err != nil {
					io.LogError(err)
					return
				}
			}
			writeAlign(al, f)
		}
		if trimMapout != "none" {
			if err = writeNameMap(namemap, trimMapout); err != nil {
				io.LogError(err)
				return
			}
		}

		if aligns.Err != nil {
			err = aligns.Err
			io.LogError(err)
		}
		return
	},
}

func writeNameMap(namemap map[string]string, outfile string) (err error) {
	var f *os.File

	if f, err = openWriteFile(outfile); err != nil {
		return
	}

	for long, short := range namemap {
		f.WriteString(long)
		f.WriteString("\t")
		f.WriteString(short)
		f.WriteString("\n")
	}
	closeWriteFile(f, outfile)
	return
}

func init() {
	trimCmd.AddCommand(nameCmd)
	nameCmd.PersistentFlags().StringVarP(&trimMapout, "out-map", "m", "none", "Mapping output file")
	nameCmd.PersistentFlags().IntVarP(&trimNb, "nb-char", "n", 1, "Number of characters to keep in sequence names")
	nameCmd.PersistentFlags().BoolVarP(&trimAuto, "auto", "a", false, "Automatically generates sequence identifiers (priority over --nb-cchar)")
}
