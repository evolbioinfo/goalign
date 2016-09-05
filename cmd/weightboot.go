package cmd

import (
	"fmt"
	"github.com/fredericlemoine/goalign/distance"
	"github.com/fredericlemoine/goalign/io"
	"github.com/spf13/cobra"
	"math/rand"
	"os"
	"time"
)

var weightbootSeed int64
var weightbootOutput string
var weightbootnb int

// weightbootCmd represents the weightboot command
var weightbootCmd = &cobra.Command{
	Use:   "weightboot",
	Short: "generate weights for all positions of the input alignment",
	Long: `generate weights for all positions of the input alignment
`,
	Run: func(cmd *cobra.Command, args []string) {
		rand.Seed(weightbootSeed)
		var f *os.File
		var err error

		if weightbootOutput == "stdout" || weightbootOutput == "-" {
			f = os.Stdout
		} else {
			f, err = os.Create(weightbootOutput)
			if err != nil {
				io.ExitWithMessage(err)
			}
		}
		for i := 0; i < weightbootnb; i++ {
			var weights []float64 = nil
			weights = distance.BuildWeights(rootalign)
			for i, w := range weights {
				if i > 0 {
					f.WriteString("\t")
				}
				f.WriteString(fmt.Sprintf("%f", w))
			}
			f.WriteString("\n")
		}
		f.Close()
	},
}

func init() {
	buildCmd.AddCommand(weightbootCmd)
	weightbootCmd.PersistentFlags().Int64VarP(&weightbootSeed, "seed", "s", time.Now().UTC().UnixNano(), "Initial Random Seed")
	weightbootCmd.PersistentFlags().StringVarP(&weightbootOutput, "output", "o", "stdout", "Weight vectors output file")
	weightbootCmd.PersistentFlags().IntVarP(&weightbootnb, "nboot", "n", 1, "Number of bootstrap replicates to build")
}
