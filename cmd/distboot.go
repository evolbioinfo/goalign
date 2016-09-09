package cmd

import (
	"fmt"
	"github.com/fredericlemoine/goalign/distance"
	"github.com/spf13/cobra"
	"math/rand"
	"os"
	"time"
)

var distbootSeed int64
var distbootOutput string
var distbootnb int
var distbootmodel string
var distbootgamma bool
var distboolRemoveGaps bool

// distbootCmd represents the distboot command
var distbootCmd = &cobra.Command{
	Use:   "distboot",
	Short: "Builds bootstrap distances matrices",
	Long: `Builds bootstrap distances matrices

If the input alignment contains several alignments, will take the first one only

Available Distances:

- pdist
- jc   : Juke-Cantor
- k2p  : Kimura 2 Parameters
- f81  : Felsenstein 81

For example:

goalign build distboot -m k2p -i align.fa -o mats.txt`,

	Run: func(cmd *cobra.Command, args []string) {
		rand.Seed(distbootSeed)
		f := openWriteFile(distbootOutput)
		align := <-rootaligns

		model := distance.Model(distbootmodel, distboolRemoveGaps)
		for i := 0; i < distbootnb; i++ {
			var weights []float64 = nil
			if distbootgamma {
				weights = distance.BuildWeights(align)
				distMatrix := distance.DistMatrix(align, weights, model, rootcpus)
				writeDistBootMatrix(distMatrix, f)
			} else {
				boot := align.BuildBootstrap()
				distMatrix := distance.DistMatrix(boot, nil, model, rootcpus)
				writeDistBootMatrix(distMatrix, f)
			}
		}
		f.Close()
	},
}

func init() {
	buildCmd.AddCommand(distbootCmd)
	distbootCmd.PersistentFlags().Int64VarP(&distbootSeed, "seed", "s", time.Now().UTC().UnixNano(), "Initial Random Seed")
	distbootCmd.PersistentFlags().StringVarP(&distbootOutput, "output", "o", "stdout", "Distance matrices output file")
	distbootCmd.PersistentFlags().StringVarP(&distbootmodel, "model", "m", "k2p", "Model for distance computation")
	distbootCmd.PersistentFlags().IntVarP(&distbootnb, "nboot", "n", 1, "Number of bootstrap replicates to build")
	distbootCmd.PersistentFlags().BoolVarP(&distbootgamma, "gamma", "g", false, "Bootstraps are done by weighting alignment using gamma generated weights")
	distbootCmd.PersistentFlags().BoolVarP(&distboolRemoveGaps, "rm-gaps", "r", false, "Do not take into account positions containing >=1 gaps")
}

func writeDistBootMatrix(matrix [][]float64, f *os.File) {
	f.WriteString(fmt.Sprintf("%d\n", len(matrix)))
	for i := 0; i < len(matrix); i++ {
		f.WriteString(fmt.Sprintf("%d", i))
		for j := 0; j < len(matrix); j++ {
			f.WriteString(fmt.Sprintf("\t%.12f", matrix[i][j]))
		}
		f.WriteString("\n")
	}
}
