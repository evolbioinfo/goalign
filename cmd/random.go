package cmd

import (
	"github.com/evolbioinfo/goalign/align"
	"github.com/evolbioinfo/goalign/io"
	"github.com/evolbioinfo/goalign/io/utils"
	"github.com/spf13/cobra"
)

var randomLength, randomSize int
var randomAA bool
var randomOutput string

// randomCmd represents the random command
var randomCmd = &cobra.Command{
	Use:   "random",
	Short: "Generate random sequences",
	Long: `Generate random sequences.
`,
	RunE: func(cmd *cobra.Command, args []string) (err error) {
		var f utils.StringWriterCloser
		var a align.Alignment

		if f, err = utils.OpenWriteFile(randomOutput); err != nil {
			io.LogError(err)
			return
		}
		defer utils.CloseWriteFile(f, randomOutput)

		if !randomAA {
			if a, err = align.RandomAlignment(align.NUCLEOTIDS, randomLength, randomSize); err != nil {
				io.LogError(err)
				return
			}
		} else {
			if a, err = align.RandomAlignment(align.AMINOACIDS, randomLength, randomSize); err != nil {
				io.LogError(err)
				return
			}
		}
		writeAlign(a, f)

		return
	},
}

func init() {
	RootCmd.AddCommand(randomCmd)
	randomCmd.PersistentFlags().IntVarP(&randomLength, "length", "l", 100, "Length of sequences to generate")
	randomCmd.PersistentFlags().IntVarP(&randomSize, "nb-seqs", "n", 10, "Number of sequences to generate")
	randomCmd.PersistentFlags().BoolVarP(&randomAA, "amino-acids", "a", false, "Aminoacid sequences (otherwise, nucleotides)")
	randomCmd.PersistentFlags().StringVarP(&randomOutput, "out-align", "o", "stdout", "Random alignment output file")
}
