package cmd

import (
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
	Run: func(cmd *cobra.Command, args []string) {
		var err error
		namemap := make(map[string]string)
		curid := 1
		f := openWriteFile(trimAlignOut)
		for al := range rootaligns.Achan {
			if rootaligns.Err != nil {
				io.ExitWithMessage(rootaligns.Err)
			}

			if trimAuto {
				if err = al.TrimNamesAuto(namemap, &curid); err != nil {
					io.ExitWithMessage(err)
				}
			} else {
				if err = al.TrimNames(namemap, trimNb); err != nil {
					io.ExitWithMessage(err)
				}
			}
			writeAlign(al, f)
		}
		writeNameMap(namemap, trimMapout)
		f.Close()
	},
}

func writeNameMap(namemap map[string]string, outfile string) {
	f := openWriteFile(outfile)
	for long, short := range namemap {
		f.WriteString(long)
		f.WriteString("\t")
		f.WriteString(short)
		f.WriteString("\n")
	}
	f.Close()
}

func init() {
	trimCmd.AddCommand(nameCmd)
	nameCmd.PersistentFlags().StringVarP(&trimMapout, "out-map", "m", "none", "Mapping output file")
	nameCmd.PersistentFlags().IntVarP(&trimNb, "nb-char", "n", 1, "Number of characters to keep in sequence names")
	nameCmd.PersistentFlags().BoolVarP(&trimAuto, "auto", "a", false, "Automatically generates sequence identifiers (priority over --nb-cchar)")
}
