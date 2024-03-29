package cmd

import (
	"fmt"

	"github.com/evolbioinfo/goalign/align"
	"github.com/evolbioinfo/goalign/io"
	"github.com/evolbioinfo/goalign/io/utils"
	"github.com/spf13/cobra"
)

var swOutput string
var swLog string
var gapopen, gapextend float64
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

Score for opening a gap is specified by --gap-open option and score for extending a gap is
specified by --gap-extend option (they should be negative).

Input file must be a fasta file containing 2 sequences. Output format may be specified
by formatting options (-p, -x, etc.)
`,
	RunE: func(cmd *cobra.Command, args []string) (err error) {
		var seqs align.SeqBag
		var al align.Alignment
		var seq1 align.Sequence
		var seq2 align.Sequence
		var ok bool
		var f, log utils.StringWriterCloser

		if f, err = utils.OpenWriteFile(swOutput); err != nil {
			io.LogError(err)
			return
		}
		defer utils.CloseWriteFile(f, swOutput)

		if swLog != "none" {
			if log, err = utils.OpenWriteFile(swLog); err != nil {
				io.LogError(err)
				return
			}
			defer utils.CloseWriteFile(log, swLog)
		}

		if seqs, err = readsequences(infile); err != nil {
			io.LogError(err)
			return
		}

		if seqs.NbSequences() != 2 {
			err = fmt.Errorf("fasta file must contain 2 sequences to align")
			io.LogError(err)
			return
		}

		if seq1, ok = seqs.Sequence(0); !ok {
			io.LogError(fmt.Errorf("sequence 0 is not present in the seqbag"))
			return
		}

		if seq2, ok = seqs.Sequence(1); !ok {
			io.LogError(fmt.Errorf("sequence 1 is not present in the seqbag"))
			return
		}

		aligner := align.NewPwAligner(seq1, seq2, align.ALIGN_ALGO_SW)
		aligner.SetGapOpenScore(gapopen)
		aligner.SetGapExtendScore(gapextend)

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
			fmt.Fprintf(log, "Align length: %d\n", aligner.Length())
			fmt.Fprintf(log, "Align Score: %.2f\n", aligner.MaxScore())
			fmt.Fprintf(log, "Align Matches: %d\n", aligner.NbMatches())
			fmt.Fprintf(log, "Align Mismatches: %d\n", aligner.NbMisMatches())
			fmt.Fprintf(log, "Align Gaps: %d\n", aligner.NbGaps())
			fmt.Fprintf(log, "Alignment:\n%s\n", aligner.AlignmentStr())
		}

		return
	},
}

func init() {
	RootCmd.AddCommand(swCmd)
	swCmd.PersistentFlags().StringVarP(&swOutput, "output", "o", "stdout", "Alignment output file")
	swCmd.PersistentFlags().StringVarP(&swLog, "log", "l", "none", "Alignment log file")
	swCmd.PersistentFlags().Float64Var(&gapopen, "gap-open", -10.0, "Score for opening a gap ")
	swCmd.PersistentFlags().Float64Var(&gapextend, "gap-extend", -0.5, "Score for extending a gap ")
	swCmd.PersistentFlags().Float64Var(&match, "match", 1.0, "Score for a match (if omitted, then take substitution matrix)")
	swCmd.PersistentFlags().Float64Var(&mismatch, "mismatch", -1.0, "Score for a mismatch (if omitted, then take substitution matrix)")
}
