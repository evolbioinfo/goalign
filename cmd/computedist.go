package cmd

import (
	"fmt"
	"log"
	"math"
	"strconv"
	"strings"

	"github.com/spf13/cobra"

	"github.com/evolbioinfo/goalign/align"
	"github.com/evolbioinfo/goalign/distance/dna"
	"github.com/evolbioinfo/goalign/distance/protein"
	"github.com/evolbioinfo/goalign/io"
	"github.com/evolbioinfo/goalign/io/utils"
	pm "github.com/evolbioinfo/goalign/models/protein"

	"gonum.org/v1/gonum/mat"
)

var computedistOutput string
var computedistModel string
var computedistRemoveGaps bool
var computedistAverage bool
var computedistAlpha float64
var computedistCountGaps int
var computedistRemoveAmbiguous bool
var computedistRange1 string
var computedistRange2 string

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

For nucleotides, differences are generally counted if nucleotides are incompatible.
For example R and Y will give a difference; N and A will not give a difference.

If distance is pdist (nucleotides), then giving the option --rm-ambiguous will not take into 
account ambiguous positions that are compatible in length normalization.
For example if --rm-ambiguous is given, then R vs. Y will be taken into account
because there is a difference. And N vs. A won't be taken into account in total length
because we are not sure whether they are identical.

For example:

goalign compute distance -m k2p -i align.ph -p
goalign compute distance -m k2p -i align.fa

if -a is given: display only the average distance

In case of a nucleotidic alignment, it is possible to specify sequence ranges to compare. For example, 
goalign compute distance -m pdist -i align.ph -p --range1 0:9 --range2 10:19
will compute distance only between sequences [0 to 9] and sequences [10 to 19].
Output matrix will be formatted the same way as usual, except that it will be made of 0 except for
the comparisons 0 vs. 10; 0 .vs 11; ...; 9 vs. 19.

`,
	RunE: func(cmd *cobra.Command, args []string) (err error) {
		var f utils.StringWriterCloser
		var model dna.DistModel
		var aligns *align.AlignChannel
		var protmodel int

		if f, err = utils.OpenWriteFile(computedistOutput); err != nil {
			io.LogError(err)
			return
		}
		defer utils.CloseWriteFile(f, computedistOutput)

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
				m.SetRemoveAmbiguous(computedistRemoveAmbiguous)
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

			var range1min, range1max, range2min, range2max int = -1, -1, -1, -1

			if computedistRange1 != "" || computedistRange2 != "" {
				r1 := strings.Split(computedistRange1, ":")
				r2 := strings.Split(computedistRange2, ":")
				if len(r1) != 2 || len(r2) != 2 {
					return fmt.Errorf("sequence ranges are not well formed, should be min:max")
				}
				if range1min, err = strconv.Atoi(r1[0]); err != nil {
					return fmt.Errorf("cannot convert range1 min to integer")
				}
				if range1max, err = strconv.Atoi(r1[1]); err != nil {
					return fmt.Errorf("cannot convert range1 max to integer")
				}
				if range2min, err = strconv.Atoi(r2[0]); err != nil {
					return fmt.Errorf("cannot convert range2 min to integer")
				}
				if range2max, err = strconv.Atoi(r2[1]); err != nil {
					return fmt.Errorf("cannot convert range2 max to integer")
				}
			}

			for align := range aligns.Achan {
				var distMatrix [][]float64
				distMatrix, err = dna.DistMatrix(align, nil, model, range1min, range1max, range2min, range2max, cmd.Flags().Changed("alpha"), computedistAlpha, rootcpus)
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
	computedistCmd.PersistentFlags().BoolVar(&computedistRemoveAmbiguous, "rm-ambiguous", false, "if true, ambiguous positions are removed for the normalisation by the length in case of non different positions. Only available for pdist (nt)")
	computedistCmd.PersistentFlags().BoolVarP(&computedistAverage, "average", "a", false, "Compute only the average distance between all pairs of sequences")
	computedistCmd.PersistentFlags().Float64Var(&computedistAlpha, "alpha", 0.0, "Gamma alpha parameter, if not given : no gamma")
	computedistCmd.PersistentFlags().StringVar(&computedistRange1, "range1", "", "If set, then will restrict distance computation to the given seq range compared to range 2 (0-based, ex --range1 0:100 means [0,100]), only for nucleotide models so far")
	computedistCmd.PersistentFlags().StringVar(&computedistRange2, "range2", "", "If set, then will restrict distance computation to the given seq range compared to range 1 (0-based, ex --range2 0:100 means [0:100]), only for nucleotide models so far")
}

func writeDistMatrix(al align.Alignment, matrix [][]float64, f utils.StringWriterCloser) (err error) {
	f.WriteString(fmt.Sprintf("%d\n", len(matrix)))
	for i := 0; i < len(matrix); i++ {
		if name, ok := al.GetSequenceNameById(i); ok {
			f.WriteString(name)
		} else {
			return fmt.Errorf("sequence %d does not exist in the alignment", i)
		}
		for j := 0; j < len(matrix); j++ {
			f.WriteString(fmt.Sprintf("\t%.12f", matrix[i][j]))
		}
		f.WriteString("\n")
	}
	return
}

func writeDistAverage(al align.Alignment, matrix [][]float64, f utils.StringWriterCloser) {
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
		log.Printf("Warning: The average distance hides %d NaN distances", nan)
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
