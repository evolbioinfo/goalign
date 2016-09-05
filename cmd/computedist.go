package cmd

import (
	"errors"
	"fmt"
	"github.com/fredericlemoine/goalign/distance"
	"github.com/fredericlemoine/goalign/io"
	"github.com/spf13/cobra"
	"os"
	"time"
)

var computedistSeed int64
var computedistOutput string
var computedistModel string

// computedistCmd represents the computedist command
var computedistCmd = &cobra.Command{
	Use:   "distance",
	Short: "Compute distance matrix of 2 sequences",
	Long: `Compute distance matrix of 2 sequences
For example:

goalign compute distance -m k2p -i align.ph -p
goalign compute distance -m k2p -i align.fa

`,
	Run: func(cmd *cobra.Command, args []string) {
		if computedistModel != "k2p" {
			io.ExitWithMessage(errors.New("Only k2p is implemented so far"))
		}
		var distMatrix [][]float64 = distance.MatrixK2P(rootalign, nil)
		writeDistMatrix(distMatrix, computedistOutput)
	},
}

func init() {
	computeCmd.AddCommand(computedistCmd)
	computedistCmd.PersistentFlags().Int64VarP(&computedistSeed, "seed", "s", time.Now().UTC().UnixNano(), "Initial Random Seed")
	computedistCmd.PersistentFlags().StringVarP(&computedistOutput, "output", "o", "stdout", "Distance matrix output file")
	computedistCmd.PersistentFlags().StringVarP(&computedistModel, "model", "m", "k2p", "Model for distance computation")
}

func writeDistMatrix(matrix [][]float64, file string) {
	var f *os.File
	var err error

	if file == "stdout" || file == "-" {
		f = os.Stdout
	} else {
		f, err = os.Create(file)
		if err != nil {
			io.ExitWithMessage(err)
		}
	}

	f.WriteString(fmt.Sprintf("%d\n", len(matrix)))
	for i := 0; i < len(matrix); i++ {
		f.WriteString(fmt.Sprintf("%d", i))
		for j := 0; j < len(matrix); j++ {
			f.WriteString(fmt.Sprintf("\t%f", matrix[i][j]))
		}
		f.WriteString("\n")
	}
	f.Close()
}
