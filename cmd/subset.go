package cmd

import (
	"bufio"
	"compress/gzip"
	"github.com/fredericlemoine/goalign/align"
	"github.com/fredericlemoine/goalign/io"
	"github.com/spf13/cobra"
	"os"
	"strings"
)

var namefile string = "stdin"
var nameout string = "stdout"

// subsetCmd represents the subset command
var subsetCmd = &cobra.Command{
	Use:   "subset",
	Short: "Take a subset of sequences from the input alignment",
	Long: `Take a subset of sequences from the input alignment

It take an alignment and a set of sequence names, and 
prints the alignments corresponding to sequence names.

For example:

goalign subset -p -i al.phy seq1 seq2 seq3
goalign subset -p -i al.phy -f seqnames_file.txt

seqnames_file should be formated with one sequence name 
per line and or coma separated. If the file contains names
 that do not exist in the alignment, they won't be taken 
into account.

The output format is the same than input format.

If -f is given, it does not take into account sequence names 
given in the comand line.
`,
	Run: func(cmd *cobra.Command, args []string) {
		subset := parseNameFile(namefile)
		out := openWriteFile(nameout)
		for al := range rootaligns {
			var filtered align.Alignment = nil
			al.Iterate(func(name string, sequence string) {
				if filtered == nil {
					filtered = align.NewAlign(align.DetectAlphabet(sequence))
				}
				_, ok := subset[name]
				if ok {
					filtered.AddSequence(name, sequence, "")
				}
			})
			writeAlign(filtered, out)
		}
		out.Close()
	},
}

func parseNameFile(file string) map[string]int {
	var f *os.File
	var r *bufio.Reader
	subset := make(map[string]int)
	var err error
	if file == "stdin" || file == "-" {
		f = os.Stdin
	} else {
		f, err = os.Open(file)
		if err != nil {
			io.ExitWithMessage(err)
		}
	}

	if strings.HasSuffix(file, ".gz") {
		if gr, err := gzip.NewReader(f); err != nil {
			io.ExitWithMessage(err)
		} else {
			r = bufio.NewReader(gr)
		}
	} else {
		r = bufio.NewReader(f)
	}

	l, e := Readln(r)
	for e == nil {
		for _, name := range strings.Split(l, ",") {
			subset[name] = 1
		}
		l, e = Readln(r)
	}
	return subset
}

func init() {
	RootCmd.AddCommand(subsetCmd)
	subsetCmd.PersistentFlags().StringVarP(&namefile, "name-file", "f", "stdin", "File containing names of sequences to keep")
	subsetCmd.PersistentFlags().StringVarP(&nameout, "output", "o", "stdout", "Alignment output file")

}
