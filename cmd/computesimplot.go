package cmd

import (
	"fmt"

	"github.com/spf13/cobra"

	"github.com/evolbioinfo/goalign/align"
	"github.com/evolbioinfo/goalign/distance/dna"
	"github.com/evolbioinfo/goalign/io"
	"github.com/evolbioinfo/goalign/io/utils"
)

var simplotOutput string
var simplotModel string
var simplotRefSeq string
var simplotWindows int

// computedistCmd represents the computedist command
var computeSimplotCmd = &cobra.Command{
	Use:   "simplot",
	Short: "Compute simplot data",
	Long: `Compute simplot data

`,
	RunE: func(cmd *cobra.Command, args []string) (err error) {
		var f utils.StringWriterCloser
		var aligns *align.AlignChannel
		var windows []struct {
			WindowStart, WindowEnd int
			CompSeq                string
			Distance               float64
		}

		if f, err = utils.OpenWriteFile(simplotOutput); err != nil {
			io.LogError(err)
			return
		}
		defer utils.CloseWriteFile(f, computedistOutput)

		if aligns, err = readalign(infile); err != nil {
			io.LogError(err)
			return
		}

		fmt.Fprintf(f, "start\tend\tcomp\tdist\n")
		for align := range aligns.Achan {

			if aligns.Err != nil {
				err = aligns.Err
				io.LogError(err)
			}
			if windows, err = dna.SimPlotDistances(align, simplotRefSeq, simplotModel, simplotWindows); err != nil {
				io.LogError(err)
				return
			}

			for _, w := range windows {
				fmt.Fprintf(f, "%d\t%d\t%s\t%v\n", w.WindowStart, w.WindowEnd, w.CompSeq, w.Distance)
			}
		}
		return
	},
}

func init() {
	computeCmd.AddCommand(computeSimplotCmd)
	computeSimplotCmd.PersistentFlags().StringVarP(&simplotOutput, "output", "o", "stdout", "Distance matrix output file")
	computeSimplotCmd.PersistentFlags().StringVarP(&simplotModel, "model", "m", "k2p", "Model for distance computation")
	computeSimplotCmd.PersistentFlags().StringVarP(&simplotRefSeq, "refseq", "r", "-", "Reference sequence to compare all others")
	computeSimplotCmd.PersistentFlags().IntVarP(&simplotWindows, "window-size", "w", 100, "Window size")
}
