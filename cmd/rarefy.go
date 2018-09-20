package cmd

import (
	"bufio"
	"compress/gzip"
	"errors"
	"os"
	"strconv"
	"strings"

	"github.com/fredericlemoine/goalign/io"
	"github.com/spf13/cobra"
)

var rarefyNb int
var rarefyOutput string
var rarefyCounts string
var rarefyReplicates int

// rarefyCmd represents the rarefy command
var rarefyCmd = &cobra.Command{
	Use:   "rarefy",
	Short: "Take a new sample taking into accounts weights",
	Long: `Take a new sample taking into accounts weights.

Each sequence in the alignment has associated counts. The sum s of the counts 
represents the number of sequences in the underlying initial dataset.

The goal is to downsample (rarefy) the initial dataset, by sampling n sequences 
from s (n<s), and taking the alignment corresponding to this new sample, i.e by 
taking only unique (different) sequences from it.

Parameters are: 
* n: the number of sequences to sample from the underlying full dataset (different
from the number of sequences in the output alignment)
* c: counts associated to each sequence (if the count of a sequence is missing, it 
is considered as 0). Sum of counts of all sequences must be > n.
* r: the number of replicates to generate (if r>1, output format will be phylip anyway)

Output: An alignment (phylip or fasta).
`,
	Run: func(cmd *cobra.Command, args []string) {
		aligns := readalign(infile)
		counts := parseCountFile(rarefyCounts)
		f := openWriteFile(rarefyOutput)
		for al := range aligns.Achan {
			if rarefyReplicates > 1 {
				rootphylip = true
			}
			for i := 0; i < rarefyReplicates; i++ {
				if sample, err := al.Rarefy(rarefyNb, counts); err != nil {
					io.ExitWithMessage(err)
				} else {
					writeAlign(sample, f)
				}
			}
		}
		f.Close()
	},
}

func init() {
	sampleCmd.AddCommand(rarefyCmd)
	rarefyCmd.PersistentFlags().IntVarP(&rarefyNb, "nb-seq", "n", 1, "Number of sequences to sample from the repeated dataset (from counts)")
	rarefyCmd.PersistentFlags().StringVarP(&rarefyOutput, "output", "o", "stdout", "Rarefied alignment output file")
	rarefyCmd.PersistentFlags().StringVarP(&rarefyCounts, "counts", "c", "stdin", "Count file (tab separated), one line per sequence: seqname\\tcount")
	rarefyCmd.PersistentFlags().IntVarP(&rarefyReplicates, "replicates", "r", 1, "Number of replicates to generate")
}

func parseCountFile(file string) map[string]int {
	var f *os.File
	var r *bufio.Reader
	counts := make(map[string]int)
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
		cols := strings.Split(l, "\t")
		if cols == nil || len(cols) != 2 {
			io.ExitWithMessage(errors.New("Bad format from counts: Wrong number of columns"))
		}
		c, err := strconv.Atoi(cols[1])
		if err != nil {
			io.ExitWithMessage(err)
		}
		counts[cols[0]] = c
		l, e = Readln(r)
	}
	return counts
}
