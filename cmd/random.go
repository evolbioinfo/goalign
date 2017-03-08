package cmd

import (
	"math/rand"
	"time"

	"github.com/fredericlemoine/goalign/align"
	"github.com/fredericlemoine/goalign/io"
	"github.com/spf13/cobra"
)

var randomLength, randomSize int
var randomAA bool
var randomOutput string
var randomSeed int64

// randomCmd represents the random command
var randomCmd = &cobra.Command{
	Use:   "random",
	Short: "Generate random sequences",
	Long: `Generate random sequences.
`,
	PersistentPreRun: func(cmd *cobra.Command, args []string) {
	},
	Run: func(cmd *cobra.Command, args []string) {
		rand.Seed(randomSeed)
		var a align.Alignment
		var err error
		f := openWriteFile(addIdOutput)
		if !randomAA {
			a, err = align.RandomAlignment(align.NUCLEOTIDS, randomLength, randomSize)
			if err != nil {
				io.ExitWithMessage(err)
			}
		} else {
			a, err = align.RandomAlignment(align.AMINOACIDS, randomLength, randomSize)
			if err != nil {
				io.ExitWithMessage(err)
			}
		}
		writeAlign(a, f)
		f.Close()
	},
}

func init() {
	RootCmd.AddCommand(randomCmd)
	randomCmd.PersistentFlags().IntVarP(&randomLength, "length", "l", 100, "Length of sequences to generate")
	randomCmd.PersistentFlags().IntVarP(&randomSize, "nb-seqs", "n", 10, "Number of sequences to generate")
	randomCmd.PersistentFlags().BoolVarP(&randomAA, "amino-acids", "a", false, "Aminoacid sequences (otherwise, nucleotides)")
	randomCmd.PersistentFlags().StringVarP(&randomOutput, "out-align", "o", "stdout", "Random alignment output file")
	randomCmd.PersistentFlags().Int64VarP(&randomSeed, "seed", "s", time.Now().UTC().UnixNano(), "Initial Random Seed")
}
