//go:build ignore

package cmd

import (
	"fmt"

	"github.com/evolbioinfo/goalign/align"
	"github.com/evolbioinfo/goalign/distance/dna"
	"github.com/evolbioinfo/goalign/io"
	"github.com/evolbioinfo/goalign/io/utils"
	"github.com/spf13/cobra"
)

var weightbootOutput string
var weightbootnb int

// weightbootCmd represents the weightboot command
var weightbootCmd = &cobra.Command{
	Use:   "weightboot",
	Short: "generate continous weights for all positions of the input alignment",
	Long: `generate continous weights for all positions of the input alignment

If the input alignment contains several alignments, will process the first one only.

Weights follow a Dirichlet distribtion D(n;1,...,1)

`,
	RunE: func(cmd *cobra.Command, args []string) (err error) {
		var f utils.StringWriterCloser
		var alignChan *align.AlignChannel

		if alignChan, err = readalign(infile); err != nil {
			io.LogError(err)
			return
		}
		al, _ := <-alignChan.Achan
		if alignChan.Err != nil {
			err = alignChan.Err
			io.LogError(err)
			return
		}

		if f, err = utils.OpenWriteFile(weightbootOutput); err != nil {
			io.LogError(err)
			return
		}
		defer utils.CloseWriteFile(f, weightbootOutput)

		for i := 0; i < weightbootnb; i++ {
			var weights []float64 = nil
			weights = dna.BuildWeightsDirichlet(al)
			for i, w := range weights {
				if i > 0 {
					f.WriteString("\t")
				}
				f.WriteString(fmt.Sprintf("%f", w))
			}
			f.WriteString("\n")
		}
		return
	},
}

func init() {
	buildCmd.AddCommand(weightbootCmd)
	weightbootCmd.PersistentFlags().StringVarP(&weightbootOutput, "output", "o", "stdout", "Weight vectors output file")
	weightbootCmd.PersistentFlags().IntVarP(&weightbootnb, "nboot", "n", 1, "Number of bootstrap replicates to build")
}
