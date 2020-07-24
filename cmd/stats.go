package cmd

import (
	"fmt"
	"os"
	"sort"

	"github.com/evolbioinfo/goalign/align"
	"github.com/evolbioinfo/goalign/io"
	"github.com/evolbioinfo/goalign/io/countprofile"
	"github.com/spf13/cobra"
)

var statpersequences bool
var statrefsequence string
var statcountprofile string

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
	If --count-profile is given along, then the output will be : unique\tnew\tboth, with:
	- 4a unique: # gaps that are unique in each sequence in the alignment
	- 4b new: # gaps that are new in each sequence compared to the profile
	- 4c both: # gaps that are unique in each sequence in the alignment and that are new compared the profile
5. Number of gap opennings (streches of gaps are counted once);
6. Number of Unique mutations;
	If --count-profile is given along, then the output will be : unique\tnew\tboth, with:
	- 6a unique: # mutations that are unique in each sequence in the alignment
	- 6b new: # mutations that are new in each sequence compared to the profile
	- 6c both: # mutations that are unique in each sequence in the alignment and that are new compared the profile
7. Number of mutations compared to a reference sequence (given with --ref-sequence, otherwise, no column);
8. Length of the sequence without gaps;
9..n Number of occurence of each character (A,C,G, etc.).

Note that --count-profile takes a tab separated file such as given by the command 
goalign stats char --per-sites

site  A C G T
0 nA  nC  nG  nT
1...
...
n...

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
				var profile *align.CountProfile
				var s string
				var ok bool
				var sb align.SeqBag
				if statcountprofile != "none" {
					if profile, err = countprofile.FromFile(statcountprofile); err != nil {
						io.LogError(err)
						return
					}
				}
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
					refseq = align.NewSequence("ref", []uint8(s), "")
				}
				err = printAllSequenceStats(al, refseq, profile)
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
	if _, ok := charmap[uint8(only[0])]; !ok && only != "*" {
		charmap[uint8(only[0])] = 0
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
		nb := charmap[uint8(k[0])]
		fmt.Fprintf(os.Stdout, "%s\t%d\t%f\n", k, nb, float64(nb)/float64(total))
	}
}

func printSiteCharStats(al align.Alignment, only string) (err error) {
	var profile *align.CountProfile
	var ok bool
	var indexonly int

	profile = align.NewCountProfileFromAlignment(al)
	onlyr := []uint8(only)
	if len(onlyr) > 1 {
		err = fmt.Errorf("Character should have length 1: %s", only)
	}

	fmt.Fprintf(os.Stdout, "site")
	if only == "*" {
		indexonly = -1
	} else {
		if indexonly, ok = profile.NameIndex(onlyr[0]); !ok {
			for site := 0; site < al.Length(); site++ {
				fmt.Fprintf(os.Stdout, "%d0%d\n", site, 0)
			}
			return
		}
	}

	for index := 0; index < profile.NbCharacters(); index++ {
		r, _ := profile.NameAt(index)
		if only == "*" || r == onlyr[0] {
			fmt.Fprintf(os.Stdout, "\t%c", r)
		}
	}
	fmt.Fprintf(os.Stdout, "\n")
	for site := 0; site < al.Length(); site++ {
		fmt.Fprintf(os.Stdout, "%d", site)
		for index := 0; index < profile.NbCharacters(); index++ {
			if indexonly == -1 || index == indexonly {
				count, _ := profile.CountAt(index, site)
				fmt.Fprintf(os.Stdout, "\t%d", count)
			}
		}
		fmt.Fprintf(os.Stdout, "\n")
	}
	return
}

func printSequenceCharStats(sb align.SeqBag, only string) (err error) {
	var sequencemap map[uint8]int

	charmap := sb.CharStats()

	// We add the only character we want to output
	// To write 0 if there are no occurences of it
	// in the alignment
	if _, ok := charmap[uint8(only[0])]; !ok && only != "*" {
		charmap[uint8(only[0])] = 0
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
				nb := sequencemap[uint8(k[0])]
				fmt.Fprintf(os.Stdout, "\t%d", nb)
			}
		}
		fmt.Fprintf(os.Stdout, "\n")
	}
	return
}

func printAllSequenceStats(al align.Alignment, refSequence align.Sequence, countProfile *align.CountProfile) (err error) {
	var sequencemap map[uint8]int

	var numnewgaps []int // new gaps that are not found in the profile
	var numnewmuts []int // new mutations that are not found in the profile

	var nummutuniques []int  // mutations that are unique in the given alignment
	var numgapsuniques []int // gaps that are unique in the given alignment

	var nummutsboth []int // mutations that are unique in the given alignment and not found in the profile
	var numgapsboth []int // gaps that are unique in the given alignment and not found in the profile

	var nummutations int
	var gaps int
	var name string
	var uniquechars []uint8

	uniquechars = al.UniqueCharacters()

	if numgapsuniques, numnewgaps, numgapsboth, err = al.NumGapsUniquePerSequence(countProfile); err != nil {
		return
	}
	if nummutuniques, numnewmuts, nummutsboth, err = al.NumMutationsUniquePerSequence(countProfile); err != nil {
		return
	}

	fmt.Fprintf(os.Stdout, "sequence")
	fmt.Fprintf(os.Stdout, "\tgaps")
	fmt.Fprintf(os.Stdout, "\tgapsstart")
	fmt.Fprintf(os.Stdout, "\tgapsend")
	fmt.Fprintf(os.Stdout, "\tgapsuniques")
	if countProfile != nil {
		fmt.Fprintf(os.Stdout, "\tgapsnew")
		fmt.Fprintf(os.Stdout, "\tgapsboth")
	}
	fmt.Fprintf(os.Stdout, "\tgapsopenning")
	fmt.Fprintf(os.Stdout, "\tmutuniques")
	if countProfile != nil {
		fmt.Fprintf(os.Stdout, "\tmutsnew")
		fmt.Fprintf(os.Stdout, "\tmutsboth")
	}
	if refSequence != nil {
		fmt.Fprintf(os.Stdout, "\tmutref")
	}
	fmt.Fprintf(os.Stdout, "\tlength")
	for _, v := range uniquechars {
		fmt.Fprintf(os.Stdout, "\t%c", v)
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
		fmt.Printf("\t%d", numgapsuniques[i])
		if countProfile != nil {
			fmt.Printf("\t%d", numnewgaps[i])
			fmt.Printf("\t%d", numgapsboth[i])
		}
		fmt.Printf("\t%d", s.NumGapsOpenning())
		fmt.Printf("\t%d", nummutuniques[i])
		if countProfile != nil {
			fmt.Printf("\t%d", numnewmuts[i])
			fmt.Printf("\t%d", nummutsboth[i])
		}
		if refSequence != nil {
			if nummutations, err = s.NumMutationsComparedToReferenceSequence(al.Alphabet(), refSequence); err != nil {
				io.LogError(err)
				return
			}
			fmt.Printf("\t%d", nummutations)
		}
		fmt.Printf("\t%d", s.Length()-gaps)
		for _, k := range uniquechars {
			nb := sequencemap[k]
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
	statsCmd.PersistentFlags().StringVar(&statcountprofile, "count-profile", "none", "A profile to compare the alignment with, and to compute statistics faster (only with --per-sequences)")
}
