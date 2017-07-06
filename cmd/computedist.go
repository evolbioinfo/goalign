package cmd

import (
	"errors"
	"fmt"
	"log"
	"math"
	"os"
	"time"

	"github.com/fredericlemoine/goalign/align"
	"github.com/fredericlemoine/goalign/distance"
	"github.com/fredericlemoine/goalign/io"
	"github.com/spf13/cobra"
)

var computedistSeed int64
var computedistOutput string
var computedistModel string
var computedistRemoveGaps bool
var computedistAverage bool

// computedistCmd represents the computedist command
var computedistCmd = &cobra.Command{
	Use:   "distance",
	Short: "Compute distance matrix of 2 sequences",
	Long: `Compute distance matrix of 2 sequences

If the input alignment contains several alignments, will compute distances
for all of them.

Available Distances:

- pdist
- jc   : Juke-Cantor
- k2p  : Kimura 2 Parameters
- f81  : Felsenstein 81
- f84  : Felsenstein 84
- tn93 : Tamura and Nei 1993

For example:

goalign compute distance -m k2p -i align.ph -p
goalign compute distance -m k2p -i align.fa

if -a is given: display only the average distance

`,
	Run: func(cmd *cobra.Command, args []string) {
		var f *os.File
		var err error

		if computedistOutput == "stdout" || computedistOutput == "-" {
			f = os.Stdout
		} else {
			f, err = os.Create(computedistOutput)
			if err != nil {
				io.ExitWithMessage(err)
			}
		}

		model := distance.Model(computedistModel, computedistRemoveGaps)
		for align := range rootaligns {
			var distMatrix [][]float64 = distance.DistMatrix(align, nil, model, rootcpus)
			if computedistAverage {
				writeDistAverage(align, distMatrix, f)
			} else {
				writeDistMatrix(align, distMatrix, f)
			}
		}
		f.Close()
	},
}

func init() {
	computeCmd.AddCommand(computedistCmd)
	computedistCmd.PersistentFlags().Int64VarP(&computedistSeed, "seed", "s", time.Now().UTC().UnixNano(), "Initial Random Seed")
	computedistCmd.PersistentFlags().StringVarP(&computedistOutput, "output", "o", "stdout", "Distance matrix output file")
	computedistCmd.PersistentFlags().StringVarP(&computedistModel, "model", "m", "k2p", "Model for distance computation")
	computedistCmd.PersistentFlags().BoolVarP(&computedistRemoveGaps, "rm-gaps", "r", false, "Do not take into account positions containing >=1 gaps")
	computedistCmd.PersistentFlags().BoolVarP(&computedistAverage, "average", "a", false, "Compute only the average distance between all pairs of sequences")
}

func writeDistMatrix(al align.Alignment, matrix [][]float64, f *os.File) {
	f.WriteString(fmt.Sprintf("%d\n", len(matrix)))
	for i := 0; i < len(matrix); i++ {
		if name, ok := al.GetSequenceNameById(i); ok {
			f.WriteString(fmt.Sprintf("%s", name))
		} else {
			io.ExitWithMessage(errors.New(fmt.Sprintf("Sequence %d does not exist in the alignment", i)))
		}
		for j := 0; j < len(matrix); j++ {
			f.WriteString(fmt.Sprintf("\t%.12f", matrix[i][j]))
		}
		f.WriteString("\n")
	}
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
