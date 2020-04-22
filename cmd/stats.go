package cmd

import (
	"fmt"
	"os"
	"sort"

	"github.com/evolbioinfo/goalign/align"
	"github.com/evolbioinfo/goalign/io"
	"github.com/spf13/cobra"
)

var statpersequences bool
var statrefsequence string

// statsCmd represents the stats command
var statsCmd = &cobra.Command{
	Use:   "stats",
	Short: "Prints different characteristics of the alignment",
	Long: `Prints different characteristics of the alignment.

1. Length of alignment;
2. Number of sequences;
3. Average number of alleles per site;
4. Number of variables sites (does ot take into account gaps or special characters);
5. Character frequencies.

If the input alignment contains several alignments, will process all of them

If --per-sequences is given, then it will print the following stats, for each sequence:

1. Number of gaps in the sequence;
2. Number of consecutive gaps at the beginning of the sequence;
3. Number of consecutive gaps at the end of the sequence;
4. Number of gaps unique to the sequence (present in no other sequence);
5. Number of gap opennings (streches of gaps are counted once);
6. Number of Unique mutations;
7. Number of mutations compared to a reference sequence (given with --ref-sequence, otherwise, no column);
8. Length of the sequence without gaps;
9..n Number of occurence of each character (A,C,G, etc.).

`,
	RunE: func(cmd *cobra.Command, args []string) (err error) {
		var aligns *align.AlignChannel

		if aligns, err = readalign(infile); err != nil {
			io.LogError(err)
			return
		}
		for al := range aligns.Achan {
			if !statpersequences {
				fmt.Fprintf(os.Stdout, "length\t%d\n", al.Length())
				fmt.Fprintf(os.Stdout, "nseqs\t%d\n", al.NbSequences())
				fmt.Fprintf(os.Stdout, "avgalleles\t%.4f\n", al.AvgAllelesPerSite())
				fmt.Fprintf(os.Stdout, "variable sites\t%d\n", al.NbVariableSites())
				printCharStats(al, "*")
				fmt.Fprintf(os.Stdout, "alphabet\t%s\n", al.AlphabetStr())
			} else {
				var refseq align.Sequence
				var s string
				var ok bool
				var sb align.SeqBag
				if statrefsequence != "none" {
					// We try to get the sequence from its name in the alignment
					if s, ok = al.GetSequence(statrefsequence); !ok {
						//Else we open the potential file
						if sb, err = readsequences(statrefsequence); err != nil {
							io.LogError(err)
							return
						}
						if sb.NbSequences() < 1 {
							err = fmt.Errorf("The reference sequence file does not contain any sequence")
							io.LogError(err)
							return
						}
						s, _ = sb.GetSequenceById(0)
					}
					refseq = align.NewSequence("ref", []rune(s), "")
				}
				err = printAllSequenceStats(al, refseq)
			}
		}

		if aligns.Err != nil {
			err = aligns.Err
			io.LogError(err)
		}
		return
	},
}

func printCharStats(align align.Alignment, only string) {
	charmap := align.CharStats()

	// We add the only character we want to output
	// To write 0 if there are no occurences of it
	// in the alignment
	if _, ok := charmap[rune(only[0])]; !ok && only != "*" {
		charmap[rune(only[0])] = 0
	}

	keys := make([]string, 0, len(charmap))
	var total int64 = 0
	for k, v := range charmap {
		if only == "*" || string(k) == only {
			keys = append(keys, string(k))
		}
		total += v
	}
	sort.Strings(keys)

	fmt.Fprintf(os.Stdout, "char\tnb\tfreq\n")
	for _, k := range keys {
		nb := charmap[rune(k[0])]
		fmt.Fprintf(os.Stdout, "%s\t%d\t%f\n", k, nb, float64(nb)/float64(total))
	}
}

func printSiteCharStats(align align.Alignment, only string) (err error) {
	var sitemap map[rune]int

	charmap := align.CharStats()

	// We add the only character we want to output
	// To write 0 if there are no occurences of it
	// in the alignment
	if _, ok := charmap[rune(only[0])]; !ok && only != "*" {
		charmap[rune(only[0])] = 0
	}

	keys := make([]string, 0, len(charmap))
	for k := range charmap {
		keys = append(keys, string(k))
	}
	sort.Strings(keys)
	fmt.Fprintf(os.Stdout, "site")
	for _, v := range keys {
		if only == "*" || v == only {
			fmt.Fprintf(os.Stdout, "\t%s", v)
		}
	}
	fmt.Fprintf(os.Stdout, "\n")
	for site := 0; site < align.Length(); site++ {
		if sitemap, err = align.CharStatsSite(site); err != nil {
			return
		}
		fmt.Fprintf(os.Stdout, "%d", site)
		for _, k := range keys {
			if only == "*" || k == only {
				nb := sitemap[rune(k[0])]
				fmt.Fprintf(os.Stdout, "\t%d", nb)
			}
		}
		fmt.Fprintf(os.Stdout, "\n")
	}
	return
}

func printSequenceCharStats(sb align.SeqBag, only string) (err error) {
	var sequencemap map[rune]int

	charmap := sb.CharStats()

	// We add the only character we want to output
	// To write 0 if there are no occurences of it
	// in the alignment
	if _, ok := charmap[rune(only[0])]; !ok && only != "*" {
		charmap[rune(only[0])] = 0
	}

	keys := make([]string, 0, len(charmap))
	for k := range charmap {
		keys = append(keys, string(k))
	}
	sort.Strings(keys)
	fmt.Fprintf(os.Stdout, "seq")
	for _, v := range keys {
		if only == "*" || v == only {
			fmt.Fprintf(os.Stdout, "\t%s", v)
		}
	}
	fmt.Fprintf(os.Stdout, "\n")
	for i := 0; i < sb.NbSequences(); i++ {
		if sequencemap, err = sb.CharStatsSeq(i); err != nil {
			return
		}
		name, _ := sb.GetSequenceNameById(i)
		fmt.Fprintf(os.Stdout, "%s", name)
		for _, k := range keys {
			if only == "*" || k == only {
				nb := sequencemap[rune(k[0])]
				fmt.Fprintf(os.Stdout, "\t%d", nb)
			}
		}
		fmt.Fprintf(os.Stdout, "\n")
	}
	return
}

func printAllSequenceStats(al align.Alignment, refSequence align.Sequence) (err error) {
	var sequencemap map[rune]int
	var numgapsunique, nummutuniques, nummutations []int
	var gaps int
	var name string

	charmap := al.CharStats()
	numgapsunique = al.NumGapsUniquePerSequence()
	nummutuniques = al.NumMutationsUniquePerSequence()
	if refSequence != nil {
		if nummutations, err = al.NumMutationsComparedToReferenceSequence(refSequence); err != nil {
			io.LogError(err)
			return
		}
	}
	keys := make([]string, 0, len(charmap))
	for k := range charmap {
		keys = append(keys, string(k))
	}
	sort.Strings(keys)
	fmt.Fprintf(os.Stdout, "sequence")
	fmt.Fprintf(os.Stdout, "\tgaps")
	fmt.Fprintf(os.Stdout, "\tgapsstart")
	fmt.Fprintf(os.Stdout, "\tgapsend")
	fmt.Fprintf(os.Stdout, "\tgapsuniques")
	fmt.Fprintf(os.Stdout, "\tgapsopenning")
	fmt.Fprintf(os.Stdout, "\tmutuniques")
	if refSequence != nil {
		fmt.Fprintf(os.Stdout, "\tmutref")
	}
	fmt.Fprintf(os.Stdout, "\tlength")
	for _, v := range keys {
		fmt.Fprintf(os.Stdout, "\t%s", v)
	}

	fmt.Fprintf(os.Stdout, "\n")
	for i, s := range al.Sequences() {
		if sequencemap, err = al.CharStatsSeq(i); err != nil {
			return
		}
		name = s.Name()
		gaps = s.NumGaps()
		fmt.Printf("%s", name)
		fmt.Printf("\t%d", gaps)
		fmt.Printf("\t%d", s.NumGapsFromStart())
		fmt.Printf("\t%d", s.NumGapsFromEnd())
		fmt.Printf("\t%d", numgapsunique[i])
		fmt.Printf("\t%d", s.NumGapsOpenning())
		fmt.Printf("\t%d", nummutuniques[i])
		if refSequence != nil {
			fmt.Printf("\t%d", nummutations[i])
		}
		fmt.Printf("\t%d", s.Length()-gaps)
		for _, k := range keys {
			nb := sequencemap[rune(k[0])]
			fmt.Printf("\t%d", nb)
		}
		fmt.Printf("\n")
	}

	return
}

// Prints the Character with the most frequency
// for each site of the alignment
func printMaxCharStats(align align.Alignment, excludeGaps bool) {
	maxchars, occur := align.MaxCharStats(excludeGaps)

	fmt.Fprintf(os.Stdout, "site\tchar\tnb\n")
	for i, c := range maxchars {
		fmt.Fprintf(os.Stdout, "%d\t%c\t%d\n", i, c, occur[i])
	}
}

func init() {
	RootCmd.AddCommand(statsCmd)
	statsCmd.PersistentFlags().BoolVar(&statpersequences, "per-sequences", false, "Prints  statistics per alignment sequences")
	statsCmd.PersistentFlags().StringVar(&statrefsequence, "ref-sequence", "none", "Reference sequence to compare each sequence with (only with --per-sequences")
}
