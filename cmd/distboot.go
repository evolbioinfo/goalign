package cmd

import (
	"fmt"

	"github.com/spf13/cobra"

	"github.com/evolbioinfo/goalign/align"
	"github.com/evolbioinfo/goalign/distance/dna"
	"github.com/evolbioinfo/goalign/distance/protein"
	"github.com/evolbioinfo/goalign/io"
	"github.com/evolbioinfo/goalign/io/utils"
	pm "github.com/evolbioinfo/goalign/models/protein"
	"gonum.org/v1/gonum/mat"
)

var distbootOutput string
var distbootnb int
var distbootAlpha float64
var distbootmodel string
var distbootcontinuous bool = false
var distbootRemoveGaps bool
var distbootFrac float64

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
Proteins:
- DAYHOFF
- JTT
- MtRev 
- LG
- WAG

For example:

goalign build distboot -m k2p -i align.fa -o mats.txt

If frac is < 1.0, bootstrap alignment used for computing distances is a partial 
bootstrap as is phylip seqboot, which means that the sites are sampled from the 
full alignment with replacement, but the bootstrap alignment length is a fraction of the
original alignment.

`,
	//If -c is given, then random continuous weights are associated to all sites.
	//Weights follow a Dirichlet distribution D(n;1,...,1)
	//`,

	RunE: func(cmd *cobra.Command, args []string) (err error) {
		var dnamodel dna.DistModel
		var protmodel *protein.ProtDistModel
		var protmodelI int
		var d *mat.Dense
		var aligns *align.AlignChannel
		var f utils.StringWriterCloser
		var weights []float64

		if aligns, err = readalign(infile); err != nil {
			io.LogError(err)
			return
		}
		if f, err = utils.OpenWriteFile(distbootOutput); err != nil {
			io.LogError(err)
			return
		}
		defer utils.CloseWriteFile(f, distbootOutput)

		align := <-aligns.Achan
		if aligns.Err != nil {
			err = aligns.Err
			io.LogError(err)
			return
		}
		if protmodelI = pm.ModelStringToInt(distbootmodel); protmodelI != -1 {
			protmodel, _ = protein.NewProtDistModel(protmodelI, true, cmd.Flags().Changed("alpha"), distbootAlpha, distbootRemoveGaps)
			protmodel.InitModel(nil, nil)
			for i := 0; i < distbootnb; i++ {
				if distbootcontinuous {
					weights = dna.BuildWeightsDirichlet(align)
					if _, _, d, err = protmodel.MLDist(align, weights); err != nil {
						io.LogError(err)
						return
					}
					writeDenseDistBootMatrix(d, align, f)
				} else {
					boot := align.BuildBootstrap(distbootFrac)
					if _, _, d, err = protmodel.MLDist(boot, nil); err != nil {
						io.LogError(err)
						return
					}
					writeDenseDistBootMatrix(d, boot, f)
				}
			}
		} else {
			if dnamodel, err = dna.Model(distbootmodel, distbootRemoveGaps); err != nil {
				io.LogError(err)
				return
			}
			for i := 0; i < distbootnb; i++ {
				var weights []float64 = nil
				var distMatrix [][]float64
				if distbootcontinuous {
					weights = dna.BuildWeightsDirichlet(align)
					if distMatrix, err = dna.DistMatrix(align, weights, dnamodel, -1, -1, -1, -1, cmd.Flags().Changed("alpha"), distbootAlpha, rootcpus); err != nil {
						io.LogError(err)
						return
					}
				} else {
					boot := align.BuildBootstrap(distbootFrac)
					if distMatrix, err = dna.DistMatrix(boot, nil, dnamodel, -1, -1, -1, -1, cmd.Flags().Changed("alpha"), distbootAlpha, rootcpus); err != nil {
						io.LogError(err)
						return
					}
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
	distbootCmd.PersistentFlags().Float64VarP(&distbootFrac, "frac", "f", 1.0, "Fraction of sites to sample (if < 1.0: Partial bootstrap as in phylip seqboot)")
	//distbootCmd.PersistentFlags().BoolVarP(&distbootcontinuous, "continuous", "c", false, "Bootstraps are done by weighting alignment with continuous weights (dirichlet)")
	distbootCmd.PersistentFlags().BoolVarP(&distbootRemoveGaps, "rm-gaps", "r", false, "Do not take into account positions containing >=1 gaps")
	distbootCmd.PersistentFlags().Float64Var(&distbootAlpha, "alpha", 0.0, "Gamma alpha parameter, if not given : no gamma")
}

func writeDistBootMatrix(matrix [][]float64, a align.Alignment, f utils.StringWriterCloser) {
	f.WriteString(fmt.Sprintf("%d\n", len(matrix)))
	for i := 0; i < len(matrix); i++ {
		name, ok := a.GetSequenceNameById(i)
		if !ok {
			f.WriteString(fmt.Sprintf("%d", i))
		} else {
			f.WriteString(name)
		}
		for j := 0; j < len(matrix); j++ {
			f.WriteString(fmt.Sprintf("\t%.12f", matrix[i][j]))
		}
		f.WriteString("\n")
	}
}

func writeDenseDistBootMatrix(matrix *mat.Dense, a align.Alignment, f utils.StringWriterCloser) {
	r, c := matrix.Dims()
	f.WriteString(fmt.Sprintf("%d\n", c))
	for i := 0; i < r; i++ {
		name, ok := a.GetSequenceNameById(i)
		if !ok {
			f.WriteString(fmt.Sprintf("%d", i))
		} else {
			f.WriteString(name)
		}
		for j := 0; j < c; j++ {
			f.WriteString(fmt.Sprintf("\t%.12f", matrix.At(i, j)))
		}
		f.WriteString("\n")
	}
}
