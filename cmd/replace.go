package cmd

import (
	"bufio"
	"compress/gzip"
	"errors"
	"fmt"
	"os"
	"strconv"
	"strings"

	"github.com/evolbioinfo/goalign/align"
	"github.com/evolbioinfo/goalign/io"
	"github.com/evolbioinfo/goalign/io/utils"
	"github.com/spf13/cobra"
)

var replaceOutput string
var replacefile string
var replaceRegexp bool
var replaceOld, replaceNew string

type repchar struct {
	seqname string
	site    int
	newchar uint8
}

// renameCmd represents the rename command
var replaceCmd = &cobra.Command{
	Use:   "replace",
	Short: "Replace characters in sequences of the input alignment (possible with a regex)",
	Long: `Replace characters in sequences of the input alignment (possible with a regex).
If the replacement changes sequence length, then returns an error.

If --posfile is given, then --old and --new are not considered. Instead, characters at sites+sequences specified in the input file 
are replaced in the alignement. The format of the input posfile is tabulated with columns:
0: sequence name
1: site index
2: new character
`,
	RunE: func(cmd *cobra.Command, args []string) (err error) {
		var aligns *align.AlignChannel
		var f utils.StringWriterCloser
		var seqs align.SeqBag
		var replace []repchar

		if f, err = utils.OpenWriteFile(replaceOutput); err != nil {
			io.LogError(err)
			return
		}
		defer utils.CloseWriteFile(f, replaceOutput)

		if unaligned {
			if !cmd.Flags().Changed("old") || !cmd.Flags().Changed("new") {
				err = errors.New("--old and --new must be specified")
				return
			}

			if seqs, err = readsequences(infile); err != nil {
				io.LogError(err)
				return
			}
			if err = seqs.Replace(replaceOld, replaceNew, replaceRegexp); err != nil {
				io.LogError(err)
				return
			}
			writeSequences(seqs, f)
		} else {
			if replacefile == "none" && (!cmd.Flags().Changed("old") || !cmd.Flags().Changed("new")) {
				err = errors.New("--old and --new must be specified")
				return
			}
			if aligns, err = readalign(infile); err != nil {
				io.LogError(err)
				return
			}

			if replacefile != "none" {
				if replace, err = readreplacefile(replacefile); err != nil {
					io.LogError(err)
					return
				}
				for al := range aligns.Achan {
					for _, rep := range replace {
						if err = al.ReplaceChar(rep.seqname, rep.site, rep.newchar); err != nil {
							io.LogError(err)
							return
						}
					}
					writeAlign(al, f)
				}
			} else {
				for al := range aligns.Achan {
					if err = al.Replace(replaceOld, replaceNew, replaceRegexp); err != nil {
						io.LogError(err)
						return
					}
					writeAlign(al, f)
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
	RootCmd.AddCommand(replaceCmd)

	replaceCmd.PersistentFlags().StringVarP(&replaceOutput, "output", "o", "stdout", "Output alignment file")
	replaceCmd.PersistentFlags().BoolVarP(&replaceRegexp, "regexp", "e", false, "Considers Replace alignment using regexp")
	replaceCmd.PersistentFlags().StringVarP(&replaceOld, "old", "s", "none", "String to replace in the sequences")
	replaceCmd.PersistentFlags().StringVarP(&replaceNew, "new", "n", "none", "New string that will replace old string in sequences")
	replaceCmd.PersistentFlags().StringVarP(&replacefile, "posfile", "f", "none", "File containing sites to replace by give characters in given sequences (deactivates --old & --new)")
	replaceCmd.PersistentFlags().BoolVar(&unaligned, "unaligned", false, "Considers input sequences as unaligned and fasta format (phylip, nexus,... options are ignored)")
}

func readreplacefile(file string) (replace []repchar, err error) {
	var f *os.File
	var r *bufio.Reader
	var gr *gzip.Reader

	var seqname string
	var site int
	var newchar uint8
	replace = make([]repchar, 0)

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
		// Authorize comments
		if !strings.HasPrefix("#") {
			cols := strings.Split(l, "\t")
			if cols == nil || len(cols) < 3 {
				err = errors.New("bad format from replace char file: There should be 3 columns: seqname\\tsite\\tnewchar")
				return
			}

			seqname = cols[0]
			if site, err = strconv.Atoi(cols[1]); err != nil {
				err = fmt.Errorf("cannot convert site index to int")
				return
			}
			newchar = uint8(cols[2][0])

			replace = append(
				replace,
				repchar{
					seqname: seqname,
					site:    site,
					newchar: newchar,
				})
		}
		l, e = utils.Readln(r)
	}
	return
}
