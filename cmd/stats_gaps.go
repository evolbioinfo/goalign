package cmd

import (
	"fmt"

	"github.com/evolbioinfo/goalign/align"
	"github.com/evolbioinfo/goalign/io"
	"github.com/evolbioinfo/goalign/io/countprofile"
	"github.com/spf13/cobra"
)

var statGapsFromStart bool
var statGapsFromEnd bool
var statGapsUnique bool
var statGapsOpenning bool
var statGapsProfile string

// charCmd represents the char command
var statGapsCmd = &cobra.Command{
	Use:   "gaps",
	Short: "Print gap stats on each alignment sequence",
	Long: `Print gap stats on each alignment sequence.

	By default, it prints, for each alignment sequence the number of gaps.

	Following options are exclusive, and given in order of priority:
	- If --from-start is specified, then counts only gaps at sequence starts;
	- If --from-end is specified, then counts only gaps at sequence ends;
	- If --unique is specified, then counts only gaps that are unique in their column. 
	  If --count-profile is given along with --unique, then the output will be : unique\tnew\tboth, with:
		- unique: # gaps that are unique in each sequence in the alignment
		- new: # gaps that are new in each sequence compared to the profile
		- both: # gaps that are unique in each sequence in the alignment and that are new compared the profile
	- If --openning is specified, then counts only gap openning (streches of gaps are counted once)
	- Otherwise, counts total number of gaps
	for the given sequence.

	Note that --count-profile takes a tab separated file such as the one given by the command 
	goalign stats char --per-sites

	site  A C G T
	0 nA  nC  nG  nT
	1...
	...
	n...
`,
	RunE: func(cmd *cobra.Command, args []string) (err error) {
		var aligns *align.AlignChannel
		var profile *align.CountProfile

		if aligns, err = readalign(infile); err != nil {
			io.LogError(err)
			return
		}

		al, _ := <-aligns.Achan
		if aligns.Err != nil {
			err = aligns.Err
			io.LogError(err)
			return
		}

		if statGapsProfile != "none" {
			if profile, err = countprofile.FromFile(statGapsProfile); err != nil {
				io.LogError(err)
				return
			}
		}

		var numnewgaps []int     // new gaps that are not found in the profile
		var numgapsuniques []int // gaps that are unique in the given alignment
		var numgapsboth []int    // gaps that are unique in the given alignment and not found in the profile

		if statGapsUnique {
			if numgapsuniques, numnewgaps, numgapsboth, err = al.NumGapsUniquePerSequence(profile); err != nil {
				io.LogError(err)
				return
			}
		}

		for i, s := range al.Sequences() {
			if statGapsFromStart {
				fmt.Printf("%s\t%d\n", s.Name(), s.NumGapsFromStart())
			} else if statGapsFromEnd {
				fmt.Printf("%s\t%d\n", s.Name(), s.NumGapsFromEnd())
			} else if statGapsUnique {
				fmt.Printf("%s\t%d", s.Name(), numgapsuniques[i])
				if statGapsProfile != "none" {
					fmt.Printf("\t%d\t%d", numnewgaps[i], numgapsboth[i])
				}
				fmt.Printf("\n")
			} else if statGapsOpenning {
				fmt.Printf("%s\t%d\n", s.Name(), s.NumGapsOpenning())
			} else {
				fmt.Printf("%s\t%d\n", s.Name(), s.NumGaps())
			}
		}

		return
	},
}

func init() {
	statGapsCmd.PersistentFlags().BoolVar(&statGapsFromStart, "from-start", false, "Count gaps in each sequence from start of sequences (until a non gap character is encountered)")
	statGapsCmd.PersistentFlags().BoolVar(&statGapsFromEnd, "from-end", false, "Count gaps in each sequence from end of sequences (until a non gap character is encountered)")
	statGapsCmd.PersistentFlags().BoolVar(&statGapsUnique, "unique", false, "Count, in each sequence, the number of gaps that are unique in a site")
	statGapsCmd.PersistentFlags().BoolVar(&statGapsOpenning, "openning", false, "Count, in each sequence, the number of gaps openning (a strech of gaps is counted once)")
	statGapsCmd.PersistentFlags().StringVar(&statGapsProfile, "count-profile", "none", "A profile to compare the alignment with, and to compute statistics faster (only with --unique)")

	statsCmd.AddCommand(statGapsCmd)
}
