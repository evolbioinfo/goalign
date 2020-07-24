package cmd

import (
	"fmt"

	"github.com/evolbioinfo/goalign/align"
	"github.com/evolbioinfo/goalign/io"
	"github.com/evolbioinfo/goalign/io/countprofile"
	"github.com/spf13/cobra"
)

var statMutationsRef string
var statMutationsUnique bool
var statMutationsProfile string

// charCmd represents the char command
var statMutationsCmd = &cobra.Command{
	Use:   "mutations",
	Short: "Print mutations stats on each alignment sequence",
	Long: `Print mutations stats on each alignment sequence.

	- If --unique is specified, then counts only mutations (characters) that are unique in their column
	for the given sequence.
	If --count-profile is given along with --unique, then the output will be : unique\tnew\tboth, with:
		- unique: # mutations that are unique in each sequence in the alignment
		- new: # mutations that are new in each sequence compared to the profile
		- both: # mutations that are unique in each sequence in the alignment and that are new compared the profile
	- If --ref-sequence is specified, it will try to extract a seqsuence having that name from the alignment. If none exist, 
	it will try to open a fasta file with the given name to take the first sequence as a reference. If a character is ambigous 
	(IUPAC notation) in an nucleotide sequence, then it is counted as a mutation only if it is incompatible with the reference character.
	

	It does not take into account '-' and 'N' as unique mutations, and does not take into account '-' and 'N' as mutations compared 
	to a reference sequence.

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

		if statMutationsProfile != "none" {
			if profile, err = countprofile.FromFile(statMutationsProfile); err != nil {
				io.LogError(err)
				return
			}
		}

		var nummutations []int
		var numnewmuts []int  // new mutations that are not found in the profile
		var nummutsboth []int // mutations that are unique in the given alignment and not found in the profile
		var num int

		if statMutationsRef != "none" {
			var s string
			var sb align.SeqBag
			var ok bool
			// We try to get the sequence from its name in the alignment
			if s, ok = al.GetSequence(statMutationsRef); !ok {
				//Else we open the potential file
				if sb, err = readsequences(statMutationsRef); err != nil {
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
			for _, s2 := range al.Sequences() {
				if num, err = s2.NumMutationsComparedToReferenceSequence(al.Alphabet(), align.NewSequence("ref", []uint8(s), "")); err != nil {
					io.LogError(err)
					return
				}
				fmt.Printf("%s\t%d\n", s2.Name(), num)
			}

		} else if statMutationsUnique {
			nummutations, numnewmuts, nummutsboth, err = al.NumMutationsUniquePerSequence(profile)
			for i, s := range al.Sequences() {
				fmt.Printf("%s\t%d", s.Name(), nummutations[i])
				if profile != nil {
					fmt.Printf("\t%d", numnewmuts[i])
					fmt.Printf("\t%d", nummutsboth[i])
				}
				fmt.Printf("\n")
			}
		} else {
			err = fmt.Errorf("Mutations should be counted by comparing to a reference sequnce with --ref-sequence")
			io.LogError(err)
			return
		}

		return
	},
}

func init() {
	statMutationsCmd.PersistentFlags().StringVar(&statMutationsRef, "ref-sequence", "none", "Reference sequence to compare each sequence with.")
	statMutationsCmd.PersistentFlags().BoolVar(&statMutationsUnique, "unique", false, "Count, in each sequence, the number of mutations/characters that are unique in a site")
	statMutationsCmd.PersistentFlags().StringVar(&statMutationsProfile, "count-profile", "none", "A profile to compare the alignment with, and to compute statistics faster (only with --unique)")

	statsCmd.AddCommand(statMutationsCmd)
}
