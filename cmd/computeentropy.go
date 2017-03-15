package cmd

import (
	"fmt"

	"github.com/spf13/cobra"
)

var entropyAverage bool

// entropyCmd represents the entropy command
var entropyCmd = &cobra.Command{
	Use:   "entropy",
	Short: "Computes entropy of a given alignment",
	Long:  `Computes entropy of a given alignment.`,
	Run: func(cmd *cobra.Command, args []string) {
		nb := 0
		if entropyAverage {
			fmt.Println("Alignment\tAvgEntropy")
		} else {
			fmt.Println("Alignment\tSite\tEntropy")
		}
		avg := 0.0
		for align := range rootaligns {
			for i := 0; i < align.Length(); i++ {
				if e, err := align.Entropy(i); err != nil {
					panic(err)
				} else {
					if entropyAverage {
						avg += e
					} else {
						fmt.Println(fmt.Sprintf("%d\t%d\t%.3f", nb, i, e))
					}
				}
			}
			if entropyAverage {
				fmt.Println(fmt.Sprintf("%d\t%.3f", nb, avg/float64(align.Length())))
			}
			nb++
		}
	},
}

func init() {
	computeCmd.AddCommand(entropyCmd)
	entropyCmd.PersistentFlags().BoolVarP(&entropyAverage, "average", "a", false, "Compute only the average entropy of input alignment")
}
