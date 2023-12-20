package cmd

import (
	"github.com/evolbioinfo/goalign/align"
	"github.com/evolbioinfo/goalign/io"
	"github.com/evolbioinfo/goalign/io/utils"
	"github.com/spf13/cobra"
)

var recombNb float64
var recombProp float64
var recombSwap bool

// recombCmd represents the recomb command
var recombCmd = &cobra.Command{
	Use:   "recomb",
	Short: "Recombine sequences in the input alignment",
	Long: `Recombine of sequences in the input alignment.

This command recombines a proportion of the input sequences into other sequences
It takes prop*nseq seqs and copy/paste a portion of them to the other prop*nseq sequences.

- if prop < 0 or prop > 0.5 : error (prop must be <= 0.5 because it will recombine x% of 
	seqs based on other x% of seqs)
- if swap is true, then swaps the two portions of sequences (2*prop sequences will be impacted)
- if swap is false, then just transfers the portion of seq1 to seq2

It may take Fasta or Phylip input alignment.

If the input alignment contains several alignments, will process all of them

Three options:
1 - The proportion of recommbining sequences. It will take n sequences 
    and will copy/paste a portion of another n sequences;
2 - The proportion of the sequence length to recombine.
3 - Swap or not

Recombine 25% of sequences by 50% (no swap):

s1 CCCCCCCCCCCCCC    s1 CCCCCCCCCCCCCC
s2 AAAAAAAAAAAAAA => s2 AAAATTTTTTTAAA
s3 GGGGGGGGGGGGGG    s3 GGGGGGGGGGGGGG
s4 TTTTTTTTTTTTTT    s4 TTTTTTTTTTTTTT

Recombine 2x25% of sequences by 50% (swap):

s1 CCCCCCCCCCCCCC    s1 CCCCCCCCCCCCCC
s2 AAAAAAAAAAAAAA => s2 AAAATTTTTTTAAA
s3 GGGGGGGGGGGGGG    s3 GGGGGGGGGGGGGG
s4 TTTTTTTTTTTTTT    s4 TTTTAAAAAAATTT

Example of usage:

goalign shuffle recomb -i align.phylip -p -n 1 -l 0.5
goalign shuffle recomb -i align.fasta -r 0.5 -n 1 -l 0.5
`,
	RunE: func(cmd *cobra.Command, args []string) (err error) {
		var aligns *align.AlignChannel
		var f utils.StringWriterCloser

		if aligns, err = readalign(infile); err != nil {
			io.LogError(err)
			return
		}
		if f, err = utils.OpenWriteFile(shuffleOutput); err != nil {
			io.LogError(err)
			return
		}
		defer utils.CloseWriteFile(f, shuffleOutput)

		for al := range aligns.Achan {
			if err = al.Recombine(recombNb, recombProp, recombSwap); err != nil {
				io.LogError(err)
				return
			}
			writeAlign(al, f)
		}

		if aligns.Err != nil {
			err = aligns.Err
			io.LogError(err)
		}
		return
	},
}

func init() {
	shuffleCmd.AddCommand(recombCmd)

	recombCmd.PersistentFlags().Float64VarP(&recombNb, "prop-seq", "n", 0.5, "Proportion of the  sequences to recombine")
	recombCmd.PersistentFlags().Float64VarP(&recombProp, "prop-length", "l", 0.5, "Proportion of length of sequences to recombine")
	recombCmd.PersistentFlags().BoolVar(&recombSwap, "swap", false, "If true, swaps sequences, otherwise just transfer seq1 subseq to seq2")
}
