package cmd

import (
	"fmt"
	"github.com/spf13/cobra"
	"os"

	"github.com/fredericlemoine/goalign/align"
	"github.com/fredericlemoine/goalign/distance"
	"github.com/fredericlemoine/goalign/io"
)

var distbootOutput string
var distbootnb int
var distbootmodel string
var distbootcontinuous bool = false
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
- f84  : Felsenstein 84
- tn93 : Tamura and Nei 1993

For example:

goalign build distboot -m k2p -i align.fa -o mats.txt
`,
	//If -c is given, then random continuous weights are associated to all sites.
	//Weights follow a Dirichlet distribution D(n;1,...,1)
	//`,

	RunE: func(cmd *cobra.Command, args []string) (err error) {
		var model distance.DistModel
		var aligns *align.AlignChannel
		var f *os.File

		if aligns, err = readalign(infile); err != nil {
			io.LogError(err)
			return
		}
		if f, err = openWriteFile(distbootOutput); err != nil {
			io.LogError(err)
			return
		}
		defer closeWriteFile(f, distbootOutput)

		align, _ := <-aligns.Achan
		if aligns.Err != nil {
			err = aligns.Err
			io.LogError(err)
			return
		}

		if model, err = distance.Model(distbootmodel, distboolRemoveGaps); err != nil {
			io.LogError(err)
			return
		}

		for i := 0; i < distbootnb; i++ {
			var weights []float64 = nil
			var distMatrix [][]float64
			if distbootcontinuous {
				weights = distance.BuildWeightsDirichlet(align)
				if distMatrix, err = distance.DistMatrix(align, weights, model, rootcpus); err != nil {
					io.LogError(err)
					return
				}
				writeDistBootMatrix(distMatrix, align, f)
			} else {
				boot := align.BuildBootstrap()
				if distMatrix, err = distance.DistMatrix(boot, nil, model, rootcpus); err != nil {
					io.LogError(err)
					return
				}
				writeDistBootMatrix(distMatrix, align, f)
			}
		}

		return
	},
}

func init() {
	buildCmd.AddCommand(distbootCmd)
	distbootCmd.PersistentFlags().StringVarP(&distbootOutput, "output", "o", "stdout", "Distance matrices output file")
	distbootCmd.PersistentFlags().StringVarP(&distbootmodel, "model", "m", "k2p", "Model for distance computation")
	distbootCmd.PersistentFlags().IntVarP(&distbootnb, "nboot", "n", 1, "Number of bootstrap replicates to build")
	//distbootCmd.PersistentFlags().BoolVarP(&distbootcontinuous, "continuous", "c", false, "Bootstraps are done by weighting alignment with continuous weights (dirichlet)")
	distbootCmd.PersistentFlags().BoolVarP(&distboolRemoveGaps, "rm-gaps", "r", false, "Do not take into account positions containing >=1 gaps")
}

func writeDistBootMatrix(matrix [][]float64, a align.Alignment, f *os.File) {
	f.WriteString(fmt.Sprintf("%d\n", len(matrix)))
	for i := 0; i < len(matrix); i++ {
		name, ok := a.GetSequenceNameById(i)
		if !ok {
			f.WriteString(fmt.Sprintf("%d", i))
		} else {
			f.WriteString(fmt.Sprintf("%s", name))
		}
		for j := 0; j < len(matrix); j++ {
			f.WriteString(fmt.Sprintf("\t%.12f", matrix[i][j]))
		}
		f.WriteString("\n")
	}
}
