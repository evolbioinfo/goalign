package cmd

import (
	"fmt"
	"github.com/fredericlemoine/goalign/io"
	"github.com/spf13/cobra"
	"math/rand"
	"os"
	"time"
)

var computedistSeed int64
var computedistOutput string
var distMatrix [][]float64

// computedistCmd represents the computedist command
var computedistCmd = &cobra.Command{
	Use:   "computedist",
	Short: "Compute distance matrix of 2 sequences",
	Long: `Compute distance matrix of 2 sequences
For example:

goalign computedist k2p -i align.ph -p
goalign computedist k2p -i align.fa

`,
	PersistentPreRun: func(cmd *cobra.Command, args []string) {
		RootCmd.PersistentPreRun(cmd, args)
		rand.Seed(computedistSeed)
	},
	PersistentPostRun: func(cmd *cobra.Command, args []string) {
		writeDistMatrix(distMatrix, computedistOutput)
	},
}

func init() {
	RootCmd.AddCommand(computedistCmd)

	computedistCmd.PersistentFlags().Int64VarP(&computedistSeed, "seed", "s", time.Now().UTC().UnixNano(), "Initial Random Seed")
	computedistCmd.PersistentFlags().StringVarP(&computedistOutput, "output", "o", "stdout", "Distance matrix output file")

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
