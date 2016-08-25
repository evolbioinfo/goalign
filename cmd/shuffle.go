package cmd

import (
	"github.com/fredericlemoine/goalign/io/fasta"
	"github.com/fredericlemoine/goalign/io/phylip"
	"github.com/spf13/cobra"
	"math/rand"
	"os"
	"time"
)

var shuffleSeed int64
var shuffleOutput string

// shuffleCmd represents the shuffle command
var shuffleCmd = &cobra.Command{
	Use:   "shuffle",
	Short: "A set of commands to shuffle an alignment",
	Long: `A set of commands to shuffle an alignment.

It takes a Fasta of Phylip alignment in input.

It is possible to:
1 - Shuffle n sites vertically: It takes n sites of the input
    alignment and reassign the characters to different sequences;
2 - Shuffle sequence order in the alignment;
3 - Recombine n sequences together.

`,
	PersistentPreRun: func(cmd *cobra.Command, args []string) {
		RootCmd.PersistentPreRun(cmd, args)
		rand.Seed(shuffleSeed)
	},
	PersistentPostRun: func(cmd *cobra.Command, args []string) {
		var f *os.File
		var err error

		if shuffleOutput == "stdout" || shuffleOutput == "-" {
			f = os.Stdout
		} else {
			f, err = os.Create(shuffleOutput)
			if err != nil {
				panic(err)
			}
		}
		if rootphylip {
			f.WriteString(phylip.WriteAlignment(rootalign))
		} else {
			f.WriteString(fasta.WriteAlignment(rootalign))
		}

		f.Close()
	},
}

func init() {
	RootCmd.AddCommand(shuffleCmd)

	shuffleCmd.PersistentFlags().Int64VarP(&shuffleSeed, "seed", "s", time.Now().UTC().UnixNano(), "Initial Random Seed")
	shuffleCmd.PersistentFlags().StringVarP(&shuffleOutput, "output", "o", "stdout", "Shuffled alignment output file")
}
