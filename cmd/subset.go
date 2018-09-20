package cmd

import (
	"bufio"
	"compress/gzip"
	"os"
	"regexp"
	"strings"

	"github.com/fredericlemoine/goalign/align"
	"github.com/fredericlemoine/goalign/io"
	"github.com/spf13/cobra"
)

var namefile string = "stdin"
var nameout string = "stdout"
var revert bool = false
var regexmatch = false

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
		var subset map[string]int
		// If input file
		if namefile != "stdin" {
			subset = parseNameFile(namefile)
		} else {
			subset = make(map[string]int)
			// Else we look at cli args
			for _, name := range args {
				subset[name] = 1
			}
		}

		out := openWriteFile(nameout)
		regexps := make([]*regexp.Regexp, 0, 10)
		if regexmatch {
			for k, _ := range subset {
				if r, err := regexp.Compile(k); err == nil {
					regexps = append(regexps, r)
				} else {
					io.ExitWithMessage(err)
				}
			}
		}

		aligns := readalign(infile)
		for al := range aligns.Achan {
			var filtered align.Alignment = nil
			al.Iterate(func(name string, sequence string) {
				if filtered == nil {
					filtered = align.NewAlign(align.DetectAlphabet(sequence))
				}
				ok := matchSeqName(name, subset, regexps, regexmatch)
				if !revert && ok {
					filtered.AddSequence(name, sequence, "")
				} else if revert && !ok {
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

// Returns true if the name is in the map
// if regexp is true then the mathc uses regexp
// otherwise, it is an exact match
func matchSeqName(name string, subset map[string]int, regexps []*regexp.Regexp, regexp bool) bool {
	ok := false
	if regexp {
		for _, r := range regexps {
			if ok = r.MatchString(name); ok {
				break
			}
		}
	} else {
		_, ok = subset[name]
	}
	return ok
}

func init() {
	RootCmd.AddCommand(subsetCmd)
	subsetCmd.PersistentFlags().StringVarP(&namefile, "name-file", "f", "stdin", "File containing names of sequences to keep")
	subsetCmd.PersistentFlags().StringVarP(&nameout, "output", "o", "stdout", "Alignment output file")
	subsetCmd.PersistentFlags().BoolVarP(&regexmatch, "regexp", "e", false, "If sequence names are given as regexp patterns")
	subsetCmd.PersistentFlags().BoolVarP(&revert, "revert", "r", false, "If true, will remove given sequences instead of keeping only them")
}
