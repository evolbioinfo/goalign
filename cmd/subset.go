package cmd

import (
	"bufio"
	"compress/gzip"
	"os"
	"regexp"
	"strconv"
	"strings"

	"github.com/evolbioinfo/goalign/align"
	"github.com/evolbioinfo/goalign/io"
	"github.com/evolbioinfo/goalign/io/utils"
	"github.com/spf13/cobra"
)

var namefile string = "stdin"
var nameout string = "stdout"
var revert bool = false
var indices bool = false
var regexmatch = false
var subsetrmgaps = false

// subsetCmd represents the subset command
var subsetCmd = &cobra.Command{
	Use:   "subset",
	Short: "Take a subset of sequences from the input alignment",
	Long: `Take a subset of sequences from the input alignment

It takes an alignment and a set of sequence names or indices, and 
prints the alignments corresponding to sequence names.

For example:

goalign subset -p -i al.phy seq1 seq2 seq3
goalign subset -p  --indices -i al.phy 0 1 2
goalign subset -p -i al.phy -f seqnames_file.txt
goalign subset -p --indices -i al.phy -f seqindices_file.txt


seqnames_file and seqindices_file should be formated with one 
sequence name/index per line and or coma separated. If the file
contains names that do not exist in the alignment, they won't be
taken into account.

The output format is the same than input format.

If -f is given, it does not take into account sequence names 
given in the comand line.
`,
	RunE: func(cmd *cobra.Command, args []string) (err error) {
		var subset map[string]int
		var aligns *align.AlignChannel
		var seqs align.SeqBag
		var f utils.StringWriterCloser
		var r *regexp.Regexp
		var indexlist []int
		var regexps []*regexp.Regexp

		// If input file
		if namefile != "stdin" {
			if subset, err = parseNameFile(namefile); err != nil {
				io.LogError(err)
				return
			}
		} else {
			subset = make(map[string]int)
			// Else we look at cli args
			for _, name := range args {
				subset[name] = 1
			}
		}

		if indices {
			var i int
			for indexstr := range subset {
				if i, err = strconv.Atoi(indexstr); err != nil {
					io.LogError(err)
					return
				}
				indexlist = append(indexlist, i)
			}
		}

		if f, err = utils.OpenWriteFile(nameout); err != nil {
			io.LogError(err)
			return
		}
		defer utils.CloseWriteFile(f, nameout)

		//regexps := make([]*regexp.Regexp, 0, 10)
		if regexmatch {
			for k := range subset {
				if r, err = regexp.Compile(k); err == nil {
					regexps = append(regexps, r)
				} else {
					io.LogError(err)
					return
				}
			}
		}

		if unaligned {
			if seqs, err = readsequences(infile); err != nil {
				io.LogError(err)
				return
			}
			var filtered align.SeqBag = nil
			var i int = 0
			seqs.Iterate(func(name string, sequence string) bool {
				if filtered == nil {
					filtered = align.NewSeqBag(seqs.Alphabet())
				}
				ok := matchSeqName(name, i, subset, regexps, regexmatch, indexlist, indices)
				if !revert && ok {
					filtered.AddSequence(name, sequence, "")
				} else if revert && !ok {
					filtered.AddSequence(name, sequence, "")
				}
				i++
				return false
			})
			writeSequences(filtered, f)
		} else {
			if aligns, err = readalign(infile); err != nil {
				io.LogError(err)
				return
			}
			for al := range aligns.Achan {
				var filtered align.Alignment = nil
				var i int = 0
				al.Iterate(func(name string, sequence string) bool {
					if filtered == nil {
						filtered = align.NewAlign(al.Alphabet())
					}
					ok := matchSeqName(name, i, subset, regexps, regexmatch, indexlist, indices)
					if !revert && ok {
						filtered.AddSequence(name, sequence, "")
					} else if revert && !ok {
						filtered.AddSequence(name, sequence, "")
					}
					i++
					return false
				})

				if subsetrmgaps {
					filtered.RemoveGapSites(1.0, false)
				}
				writeAlign(filtered, f)
			}
			if aligns.Err != nil {
				err = aligns.Err
				io.LogError(err)
			}
		}
		return
	},
}

func parseNameFile(file string) (subset map[string]int, err error) {
	var f *os.File
	var r *bufio.Reader
	var gr *gzip.Reader

	subset = make(map[string]int)

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
		} else {
			r = bufio.NewReader(gr)
		}
	} else {
		r = bufio.NewReader(f)
	}

	l, e := utils.Readln(r)
	for e == nil {
		for _, name := range strings.Split(l, ",") {
			subset[name] = 1
		}
		l, e = utils.Readln(r)
	}
	return
}

// Returns true if the name is in the map
// if regexp is true then the mathc uses regexp
// otherwise, it is an exact match
func matchSeqName(name string, index int, subset map[string]int, regexps []*regexp.Regexp, regexp bool, indexlist []int, indices bool) bool {
	ok := false
	if regexp {
		for _, r := range regexps {
			if ok = r.MatchString(name); ok {
				break
			}
		}
	} else if indices {
		for _, r := range indexlist {
			if ok = (r == index); ok {
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
	subsetCmd.PersistentFlags().BoolVarP(&regexmatch, "regexp", "e", false, "If sequence names are given as regexp patterns (has priority over --indices)")
	subsetCmd.PersistentFlags().BoolVarP(&revert, "revert", "r", false, "If true, will remove given sequences instead of keeping only them")
	subsetCmd.PersistentFlags().BoolVar(&indices, "indices", false, "If true, extracts given sequence indices instead of sequence names (0-based)")
	subsetCmd.PersistentFlags().BoolVar(&unaligned, "unaligned", false, "Considers input sequences as unaligned and fasta format (phylip, nexus,... options are ignored)")
	subsetCmd.PersistentFlags().BoolVar(&subsetrmgaps, "remove-gaps", false, "If true, then remove gap only sites of the sub-alignment (only when --unaligned is not given)")
}
