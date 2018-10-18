package cmd

import (
	"fmt"
	"os"

	"github.com/fredericlemoine/goalign/align"
	"github.com/fredericlemoine/goalign/io"
	"github.com/spf13/cobra"
)

var swOutput string
var swLog string
var gap float64
var match float64
var mismatch float64

// translateCmd represents the addid command
var swCmd = &cobra.Command{
	Use:   "sw",
	Short: "Aligns 2 sequences using Smith&Waterman algorithm",
	Long: `Aligns 2 sequences using Smith&Waterman algorithm.

Input : Fasta file
Output: Aligned file (format depending on format options)

If neither --match nor --mismatch are specified, then match and mismatch scores
are taken from blosum62 or dnafull substitution matrices (taken from EMBOSS WATER)
depending on the input sequences alphabets.

Only one kind of gap penalty is considered so far (no gap extension).
`,
	RunE: func(cmd *cobra.Command, args []string) (err error) {
		var seqs align.SeqBag
		var al align.Alignment
		var seq1 align.Sequence
		var seq2 align.Sequence
		var ok bool
		var f, log *os.File

		if f, err = openWriteFile(swOutput); err != nil {
			io.LogError(err)
			return
		}
		defer closeWriteFile(f, swOutput)

		if swLog != "none" {
			if log, err = openWriteFile(swLog); err != nil {
				io.LogError(err)
				return
			}
			defer closeWriteFile(log, swLog)
		}

		if seqs, err = readsequences(infile); err != nil {
			io.LogError(err)
			return
		}

		if seqs.NbSequences() != 2 {
			err = fmt.Errorf("Fasta file must contain 2 sequences to align")
			io.LogError(err)
			return
		}

		if seq1, ok = seqs.Sequence(0); !ok {
			io.LogError(fmt.Errorf("Sequence 0 is not present in the seqbag"))
			return
		}

		if seq2, ok = seqs.Sequence(1); !ok {
			io.LogError(fmt.Errorf("Sequence 1 is not present in the seqbag"))
			return
		}

		aligner := align.NewPwAligner(seq1, seq2, align.ALIGN_ALGO_SW)
		aligner.SetGapScore(gap)

		if cmd.Flags().Changed("mismatch") || cmd.Flags().Changed("match") {
			aligner.SetScore(match, mismatch)
		}

		if al, err = aligner.Alignment(); err != nil {
			io.LogError(err)
			return
		}
		writeAlign(al, f)

		if log != nil {
			start1, start2 := aligner.AlignStarts()
			end1, end2 := aligner.AlignEnds()
			fmt.Fprintf(log, "Query Start,End: %d,%d\n", start1, end1)
			fmt.Fprintf(log, "Subject Start,End: %d,%d\n", start2, end2)
			fmt.Fprintf(log, "Align length: %d\n", al.Length())
		}

		return
	},
}

func init() {
	RootCmd.AddCommand(swCmd)
	swCmd.PersistentFlags().StringVarP(&swOutput, "output", "o", "stdout", "Alignment output file")
	swCmd.PersistentFlags().StringVarP(&swLog, "log", "l", "none", "Alignment log file")
	swCmd.PersistentFlags().Float64Var(&gap, "gap", -1.0, "Score for a gap")
	swCmd.PersistentFlags().Float64Var(&match, "match", 1.0, "Score for a match (if omitted, then take substitution matrix)")
	swCmd.PersistentFlags().Float64Var(&mismatch, "mismatch", -1.0, "Score for a mismatch (if omitted, then take substitution matrix)")
}
