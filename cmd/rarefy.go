package cmd

import (
	"bufio"
	"compress/gzip"
	"errors"
	"os"
	"strconv"
	"strings"

	"github.com/evolbioinfo/goalign/align"
	"github.com/evolbioinfo/goalign/io"
	"github.com/evolbioinfo/goalign/io/utils"
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
	RunE: func(cmd *cobra.Command, args []string) (err error) {
		var f *os.File
		var counts map[string]int

		if f, err = openWriteFile(rarefyOutput); err != nil {
			io.LogError(err)
			return
		}
		defer closeWriteFile(f, rarefyOutput)

		if counts, err = parseCountFile(rarefyCounts); err != nil {
			io.LogError(err)
			return
		}

		if unaligned {
			var seqs align.SeqBag
			var sample align.SeqBag

			if seqs, err = readsequences(infile); err != nil {
				io.LogError(err)
				return
			}
			if rarefyReplicates > 1 {
				rootphylip = true
			}
			for i := 0; i < rarefyReplicates; i++ {
				if sample, err = seqs.RarefySeqBag(rarefyNb, counts); err != nil {
					io.LogError(err)
					return
				}
				writeSequences(sample, f)
			}
		} else {
			var aligns *align.AlignChannel
			var sample align.Alignment

			if aligns, err = readalign(infile); err != nil {
				io.LogError(err)
				return
			}

			for al := range aligns.Achan {
				if rarefyReplicates > 1 {
					rootphylip = true
				}
				for i := 0; i < rarefyReplicates; i++ {
					if sample, err = al.Rarefy(rarefyNb, counts); err != nil {
						io.LogError(err)
						return
					}
					writeAlign(sample, f)
				}
			}

			if aligns.Err != nil {
				err = aligns.Err
				io.LogError(err)
			}
		}

		return
	},
}

func init() {
	sampleCmd.AddCommand(rarefyCmd)
	rarefyCmd.PersistentFlags().BoolVar(&unaligned, "unaligned", false, "Considers sequences as unaligned and format fasta (phylip, nexus,... options are ignored)")
	rarefyCmd.PersistentFlags().IntVarP(&rarefyNb, "nb-seq", "n", 1, "Number of sequences to sample from the repeated dataset (from counts)")
	rarefyCmd.PersistentFlags().StringVarP(&rarefyOutput, "output", "o", "stdout", "Rarefied alignment output file")
	rarefyCmd.PersistentFlags().StringVarP(&rarefyCounts, "counts", "c", "stdin", "Count file (tab separated), one line per sequence: seqname\\tcount")
	rarefyCmd.PersistentFlags().IntVarP(&rarefyReplicates, "replicates", "r", 1, "Number of replicates to generate")
}

func parseCountFile(file string) (counts map[string]int, err error) {
	var f *os.File
	var r *bufio.Reader
	var gr *gzip.Reader
	var c int

	counts = make(map[string]int)

	if file == "stdin" || file == "-" {
		f = os.Stdin
	} else {
		if f, err = os.Open(file); err != nil {
			return
		}
	}

	if strings.HasSuffix(file, ".gz") {
		if gr, err = gzip.NewReader(f); err != nil {
			return
		}
		r = bufio.NewReader(gr)
	} else {
		r = bufio.NewReader(f)
	}
	l, e := utils.Readln(r)
	for e == nil {
		cols := strings.Split(l, "\t")
		if cols == nil || len(cols) != 2 {
			err = errors.New("Bad format from counts: Wrong number of columns")
			return
		}
		if c, err = strconv.Atoi(cols[1]); err != nil {
			return
		}
		counts[cols[0]] = c
		l, e = utils.Readln(r)
	}

	return
}
