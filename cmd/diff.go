package cmd

import (
	"fmt"
	"os"
	"sort"
	"strings"

	"github.com/evolbioinfo/goalign/align"
	"github.com/evolbioinfo/goalign/io"
	"github.com/spf13/cobra"
)

var diffOutput string
var diffCount bool
var diffNoGaps bool
var reverse bool

// statsCmd represents the stats command
var diffCmd = &cobra.Command{
	Use:   "diff",
	Short: "Prints only characters that are different from the first sequence of the alignment",
	Long: `Prints only characters that are different from the first sequence of the alignment.

Takes an input alignment, and compares all the sequences to the first one.
Any character that is identical to the reference sequence is replaced with ".".

If option --counts is given, then the output is not an alignment but a count file containing,
for each sequence, the number of occurence of each difference with the reference sequence. 
The format is tab separated, with following columns:

1. Sequence name (reference sequence is not included)
2,...,end: For each type of change, its number of occurence


If option --reverse is given, then replaces . with the characters on the first sequence
`,
	RunE: func(cmd *cobra.Command, args []string) (err error) {
		var aligns *align.AlignChannel
		var f *os.File

		if f, err = openWriteFile(diffOutput); err != nil {
			io.LogError(err)
			return
		}
		defer closeWriteFile(f, diffOutput)

		if aligns, err = readalign(infile); err != nil {
			io.LogError(err)
			return
		}

		for al := range aligns.Achan {
			if diffCount {
				alldiffs, diffs := al.CountDifferences()
				writeDiffCounts(al, alldiffs, diffs, f)
			} else {
				if reverse {
					al.ReplaceMatchChars()
				} else {
					al.DiffWithFirst()
				}
				writeAlign(al, f)
			}
		}

		if aligns.Err != nil {
			err = aligns.Err
			io.LogError(err)
		}
		return
	},
}

func writeDiffCounts(al align.Alignment, alldiffs []string, diffs []map[string]int, f *os.File) {
	sort.Strings(alldiffs)
	for _, d := range alldiffs {
		if !(strings.Contains(d, "-") && diffNoGaps) {
			fmt.Fprintf(f, "\t%s", d)
		}
	}
	fmt.Fprintf(f, "\n")

	for i := 0; i < len(diffs); i++ {
		name, _ := al.GetSequenceNameById(i + 1)
		fmt.Fprintf(f, "%s", name)
		for _, d := range alldiffs {
			if !(strings.Contains(d, "-") && diffNoGaps) {
				fmt.Fprintf(f, "\t%d", diffs[i][d])
			}
		}
		fmt.Fprintf(f, "\n")
	}
}

func init() {
	diffCmd.PersistentFlags().StringVarP(&diffOutput, "output", "o", "stdout", "Diff output file")
	diffCmd.PersistentFlags().BoolVar(&diffCount, "counts", false, "Count differences instead of writting only identical characters")
	diffCmd.PersistentFlags().BoolVar(&diffNoGaps, "no-gaps", false, "Do not count gaps (only with --counts)")
	diffCmd.PersistentFlags().BoolVar(&reverse, "reverse", false, "Restore identical characters (.) using first sequence positions")
	RootCmd.AddCommand(diffCmd)
}
