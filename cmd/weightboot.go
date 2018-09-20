// +build ignore

package cmd

import (
	"fmt"
	"os"

	"github.com/fredericlemoine/goalign/align"
	"github.com/fredericlemoine/goalign/distance"
	"github.com/fredericlemoine/goalign/io"
	"github.com/spf13/cobra"
)

// weightbootCmd represents the weightboot command
var weightbootCmd = &cobra.Command{
	Use:   "weightboot",
	Short: "generate continous weights for all positions of the input alignment",
	Long: `generate continous weights for all positions of the input alignment

If the input alignment contains several alignments, will process the first one only.

Weights follow a Dirichlet distribtion D(n;1,...,1)

`,
	RunE: func(cmd *cobra.Command, args []string) (err error) {
		var f *os.File
		var aligns align.AlignChannel

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

		if f, err = openWriteFile(weightbootOutput); err != nil {
			io.LogError(err)
			return
		}
		defer closeWriteFile(f, weightbootOutput)

		for i := 0; i < weightbootnb; i++ {
			var weights []float64 = nil
			weights = distance.BuildWeightsDirichlet(al)
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
