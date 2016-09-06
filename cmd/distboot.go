package cmd

import (
	"errors"
	"fmt"
	"github.com/fredericlemoine/goalign/distance"
	"github.com/fredericlemoine/goalign/io"
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

// distbootCmd represents the distboot command
var distbootCmd = &cobra.Command{
	Use:   "distboot",
	Short: "Builds bootstrap distances matrices",
	Long: `Builds bootstrap distances matrices

If the input alignment contains several alignments, will take the first one only

For example:

goalign build distboot -m k2p -i align.fa -o mats.txt`,

	Run: func(cmd *cobra.Command, args []string) {
		rand.Seed(distbootSeed)
		if distbootmodel != "k2p" {
			io.ExitWithMessage(errors.New("Only k2p is implemented so far"))
		}
		var f *os.File
		var err error

		if distbootOutput == "stdout" || distbootOutput == "-" {
			f = os.Stdout
		} else {
			f, err = os.Create(distbootOutput)
			if err != nil {
				io.ExitWithMessage(err)
			}
		}

		align := <-rootaligns

		for i := 0; i < distbootnb; i++ {
			var weights []float64 = nil
			if distbootgamma {
				weights = distance.BuildWeights(align)
				distMatrix := distance.MatrixK2P(align, weights)
				writeDistBootMatrix(distMatrix, f)
			} else {
				boot := align.BuildBootstrap()
				distMatrix := distance.MatrixK2P(boot, nil)
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

}

func writeDistBootMatrix(matrix [][]float64, f *os.File) {
	f.WriteString(fmt.Sprintf("%d\n", len(matrix)))
	for i := 0; i < len(matrix); i++ {
		f.WriteString(fmt.Sprintf("%d", i))
		for j := 0; j < len(matrix); j++ {
			f.WriteString(fmt.Sprintf("\t%f", matrix[i][j]))
		}
		f.WriteString("\n")
	}
}
