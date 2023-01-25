package cmd

import (
	"fmt"

	"github.com/evolbioinfo/goalign/align"
	"github.com/evolbioinfo/goalign/io"
	"github.com/spf13/cobra"
)

var statMutationsListAA bool

// charCmd represents the char command
var statMutationsListCmd = &cobra.Command{
	Use:   "list",
	Short: "Print mutation list of each alignment sequence",
	Long: `Print mutations list of each alignment sequence.

	- --ref-sequence: it will try to extract a sequence having that name from the alignment. If none exist, 
	it will try to open a fasta file with the given name to take the first sequence as a reference. If a character is ambigous 
	(IUPAC notation) in an nucleotide sequence, then it is counted as a mutation only if it is incompatible with the reference character.

	It does not take into account 'N' as mutations compared to a reference sequence.
`,
	RunE: func(cmd *cobra.Command, args []string) (err error) {
		var aligns *align.AlignChannel
		var mutations []align.Mutation

		if aligns, err = readalign(infile); err != nil {
			io.LogError(err)
			return
		}

		al := <-aligns.Achan
		if aligns.Err != nil {
			err = aligns.Err
			io.LogError(err)
			return
		}

		if statMutationsRef == "none" {
			err = fmt.Errorf("mutations should be counted by comparing to a reference sequnce with --ref-sequence")
			io.LogError(err)
			return
		}

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
				err = fmt.Errorf("the reference sequence file does not contain any sequence")
				io.LogError(err)
				return
			}
			s, _ = sb.GetSequenceById(0)
		}
		for _, s2 := range al.Sequences() {
			if s2.Name() != statMutationsRef {
				if mutations, err = s2.ListMutationsComparedToReferenceSequence(al.Alphabet(), align.NewSequence("ref", []uint8(s), ""), statMutationsListAA); err != nil {
					io.LogError(err)
					return
				}
				fmt.Printf("%s", s2.Name())
				for i, m := range mutations {
					if i == 0 {
						fmt.Printf("\t%c%d%s", m.Ref, m.Pos, string(m.Alt))
					} else {
						fmt.Printf(",%c%d%s", m.Ref, m.Pos, string(m.Alt))
					}
				}
				fmt.Printf("\n")
			}
		}
		return
	},
}

func init() {
	statMutationsListCmd.PersistentFlags().BoolVar(&statMutationsListAA, "aa", false, "Take the reference sequence condon by codon, and translate ")
	statMutationsCmd.AddCommand(statMutationsListCmd)
}
