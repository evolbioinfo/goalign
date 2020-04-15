package cmd

import (
	"fmt"

	"github.com/evolbioinfo/goalign/align"
	"github.com/evolbioinfo/goalign/io"
	"github.com/spf13/cobra"
)

var statMutationsRef string
var statMutationsUnique bool

// charCmd represents the char command
var statMutationsCmd = &cobra.Command{
	Use:   "mutations",
	Short: "Print mutations stats on each alignment sequence",
	Long: `Print mutations stats on each alignment sequence.

	- If --unique is specified, then counts only mutations (characters) that are unique in their column
	for the given sequence.
	- If --ref-sequence is specified, it will try to extract a seqsuence having that name from the alignment. If none exist, 
	it will try to open a fasta file with the given name to take the first sequence as a reference. If a character is ambigous 
	(IUPAC notation) in an nucleotide sequence, then it is counted as a mutation only if it is incompatible with the reference character.
	

	It does not take into account '-' and 'N' as unique mutations, and does not take into account '-' and 'N' as mutations compared 
	to a reference sequence.
`,
	RunE: func(cmd *cobra.Command, args []string) (err error) {
		var aligns *align.AlignChannel

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

		var nummutations []int
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
			if nummutations, err = al.NumMutationsComparedToReferenceSequence(align.NewSequence("ref", []rune(s), "")); err != nil {
				io.LogError(err)
				return
			}
		} else if statMutationsUnique {
			nummutations = al.NumMutationsUniquePerSequence()
		} else {
			err = fmt.Errorf("Mutations should be counted by comparing to a reference sequnce with --ref-sequence")
			io.LogError(err)
			return
		}
		for i, s := range al.Sequences() {
			fmt.Printf("%s\t%d\n", s.Name(), nummutations[i])
		}

		return
	},
}

func init() {
	statMutationsCmd.PersistentFlags().StringVar(&statMutationsRef, "ref-sequence", "none", "Reference sequence to compare each sequence with.")
	statMutationsCmd.PersistentFlags().BoolVar(&statMutationsUnique, "unique", false, "Count, in each sequence, the number of mutations/characters that are unique in a site")

	statsCmd.AddCommand(statMutationsCmd)
}
