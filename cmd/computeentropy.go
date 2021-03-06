package cmd

import (
	"fmt"
	"math"

	"github.com/evolbioinfo/goalign/align"
	"github.com/evolbioinfo/goalign/io"
	"github.com/spf13/cobra"
)

var entropyAverage bool
var entropyRemoveGaps bool

// entropyCmd represents the entropy command
var entropyCmd = &cobra.Command{
	Use:   "entropy",
	Short: "Computes entropy of a given alignment",
	Long: `Computes entropy of a given alignment.

Example: 
goalign compute entropy -i alignment.fa
goalign compute entropy -i alignment.phy -p

It is possible to compute the average entropy:
goalign compute entropy -i alignment.phy -p -a

Which will print one average entropy per alignment in the input file:
Alignment \t AvgEntropy

Otherwise, it will print one entropy per alignment site, in a tab separated form:
Alignment \t Site \t Entropy

the computation does not take into account the following characters:
-> '*'
-> '-' (if --remove-gaps is given)

If a site is made fully of '-' (if --remove-gaps is given) or '*', then its entropy will be "NaN",
and it will not be taken into account in the average.
`,
	RunE: func(cmd *cobra.Command, args []string) (err error) {
		var aligns *align.AlignChannel
		var e float64

		if aligns, err = readalign(infile); err != nil {
			io.LogError(err)
			return
		}

		nb := 0
		if entropyAverage {
			fmt.Println("Alignment\tAvgEntropy")
		} else {
			fmt.Println("Alignment\tSite\tEntropy")
		}
		for align := range aligns.Achan {
			avg := 0.0
			total := 0
			for i := 0; i < align.Length(); i++ {
				if e, err = align.Entropy(i, entropyRemoveGaps); err != nil {
					io.LogError(err)
					return
				} else {
					if entropyAverage {
						if !math.IsNaN(e) {
							avg += e
							total++
						}
					} else {
						fmt.Printf("%d\t%d\t%.3f\n", nb, i, e)
					}
				}
			}
			if entropyAverage {
				fmt.Printf("%d\t%.3f\n", nb, avg/float64(total))
			}
			nb++
		}

		if aligns.Err != nil {
			err = aligns.Err
			io.LogError(err)
		}
		return
	},
}

func init() {
	computeCmd.AddCommand(entropyCmd)
	entropyCmd.PersistentFlags().BoolVarP(&entropyAverage, "average", "a", false, "Compute only the average entropy of input alignment")
	entropyCmd.PersistentFlags().BoolVarP(&entropyRemoveGaps, "remove-gaps", "g", false, "If true, then do not take into account gaps in the computation")
}
