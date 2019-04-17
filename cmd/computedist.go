package cmd

import (
	"errors"
	"fmt"
	"github.com/spf13/cobra"
	"log"
	"math"
	"os"

	"github.com/evolbioinfo/goalign/align"
	"github.com/evolbioinfo/goalign/distance/dna"
	"github.com/evolbioinfo/goalign/distance/protein"
	"github.com/evolbioinfo/goalign/io"
	pm "github.com/evolbioinfo/gotree/models/protein"

	"gonum.org/v1/gonum/mat"
)

var computedistOutput string
var computedistModel string
var computedistRemoveGaps bool
var computedistAverage bool
var computedistAlpha float64
var computedistCountGaps int

// computedistCmd represents the computedist command
var computedistCmd = &cobra.Command{
	Use:   "distance",
	Short: "Compute distance matrix ",
	Long: `Compute distance matrix

If the input alignment contains several alignments, will compute distances
for all of them.

Available Distances:

Nucleotides:
- pdist
- rawdist : raw distance (like pdist, without normalization by length)
- jc      : Juke-Cantor
- k2p     : Kimura 2 Parameters
- f81     : Felsenstein 81
- f84     : Felsenstein 84
- tn93    : Tamura and Nei 1993
Proteins:
- DAYHOFF
- JTT
- MtRev 
- LG
- WAG

For example:

goalign compute distance -m k2p -i align.ph -p
goalign compute distance -m k2p -i align.fa

if -a is given: display only the average distance

`,
	RunE: func(cmd *cobra.Command, args []string) (err error) {
		var f *os.File
		var model dna.DistModel
		var aligns *align.AlignChannel
		var protmodel int

		if f, err = openWriteFile(computedistOutput); err != nil {
			io.LogError(err)
			return
		}
		defer closeWriteFile(f, computedistOutput)

		if aligns, err = readalign(infile); err != nil {
			io.LogError(err)
			return
		}

		// If prot model
		if protmodel = pm.ModelStringToInt(computedistModel); protmodel != -1 {
			var d *mat.Dense
			m, _ := protein.NewProtDistModel(protmodel, true, cmd.Flags().Changed("alpha"), computedistAlpha, computedistRemoveGaps)
			m.InitModel(nil, nil)
			for align := range aligns.Achan {
				if _, _, d, err = m.MLDist(align, nil); err != nil {
					io.LogError(err)
					return
				}
				if computedistAverage {
					writeDistAverage(align, denseToSlice(d), f)
				} else {
					if err = writeDistMatrix(align, denseToSlice(d), f); err != nil {
						io.LogError(err)
						return
					}
				}
			}

			if aligns.Err != nil {
				err = aligns.Err
				io.LogError(err)
			}

		} else {
			switch computedistModel {
			case "rawdist":
				m := dna.NewRawDistModel(computedistRemoveGaps)
				m.SetCountGapMutations(computedistCountGaps)
				if err = m.SetCountGapMutations(computedistCountGaps); err != nil {
					io.LogError(err)
					return
				}
				model = m
			case "pdist":
				m := dna.NewPDistModel(computedistRemoveGaps)
				if err = m.SetCountGapMutations(computedistCountGaps); err != nil {
					io.LogError(err)
					return
				}
				model = m
			default:
				if model, err = dna.Model(computedistModel, computedistRemoveGaps); err != nil {
					io.LogError(err)
					return
				}
			}

			for align := range aligns.Achan {
				var distMatrix [][]float64
				distMatrix, err = dna.DistMatrix(align, nil, model, cmd.Flags().Changed("alpha"), computedistAlpha, rootcpus)
				if err != nil {
					io.LogError(err)
					return
				}
				if computedistAverage {
					writeDistAverage(align, distMatrix, f)
				} else {
					if err = writeDistMatrix(align, distMatrix, f); err != nil {
						io.LogError(err)
						return
					}
				}
			}

			if aligns.Err != nil {
				err = aligns.Err
				io.LogError(err)
			}
		}
		return
	},
}

func init() {
	computeCmd.AddCommand(computedistCmd)
	computedistCmd.PersistentFlags().StringVarP(&computedistOutput, "output", "o", "stdout", "Distance matrix output file")
	computedistCmd.PersistentFlags().StringVarP(&computedistModel, "model", "m", "k2p", "Model for distance computation")
	computedistCmd.PersistentFlags().BoolVarP(&computedistRemoveGaps, "rm-gaps", "r", false, "Do not take into account positions containing >=1 gaps")
	computedistCmd.PersistentFlags().IntVar(&computedistCountGaps, "gap-mut", 0, "Count gaps to nt as mutations: 0: inactivated, 1: only internal gaps, 2: all gaps. Only available for rawdist and pdist (nt)")
	computedistCmd.PersistentFlags().BoolVarP(&computedistAverage, "average", "a", false, "Compute only the average distance between all pairs of sequences")
	computedistCmd.PersistentFlags().Float64Var(&computedistAlpha, "alpha", 0.0, "Gamma alpha parameter, if not given : no gamma")
}

func writeDistMatrix(al align.Alignment, matrix [][]float64, f *os.File) (err error) {
	f.WriteString(fmt.Sprintf("%d\n", len(matrix)))
	for i := 0; i < len(matrix); i++ {
		if name, ok := al.GetSequenceNameById(i); ok {
			f.WriteString(fmt.Sprintf("%s", name))
		} else {
			return errors.New(fmt.Sprintf("Sequence %d does not exist in the alignment", i))
		}
		for j := 0; j < len(matrix); j++ {
			f.WriteString(fmt.Sprintf("\t%.12f", matrix[i][j]))
		}
		f.WriteString("\n")
	}
	return
}

func writeDistAverage(al align.Alignment, matrix [][]float64, f *os.File) {
	sum := 0.0
	total := 0
	nan := 0
	for i := 0; i < len(matrix); i++ {
		for j := i + 1; j < len(matrix); j++ {
			if math.IsNaN(matrix[i][j]) {
				nan++
			} else {
				sum += matrix[i][j]
				total++
			}
		}
	}
	if nan > 0 {
		log.Print(fmt.Sprintf("Warning: The average distance hides %d NaN distances", nan))
	}
	sum /= float64(total)
	f.WriteString(fmt.Sprintf("%.12f\n", sum))
}

func denseToSlice(m *mat.Dense) (s [][]float64) {
	r, c := m.Dims()
	s = make([][]float64, r)
	for i := 0; i < r; i++ {
		s[i] = make([]float64, c)
		for j := 0; j < c; j++ {
			s[i][j] = m.At(i, j)
		}
	}
	return
}
