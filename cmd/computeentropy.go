package cmd

import (
	"fmt"
	"math"
	"strconv"

	"github.com/evolbioinfo/goalign/align"
	"github.com/evolbioinfo/goalign/io"
	"github.com/spf13/cobra"
)

var entropyAverage bool
var entropyRemoveGaps bool
var entropyRemoveAmbiguities bool
var entropyLogBase string
var entropyPrecision int

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

Entropy is Shannon entropy, computed with the natural logarithm by default, i.e. values
are in nats. Use --log-base 2 to get values in bits (as in "goalign compute pssm -n 4"
sequence logos, which are always in bits), or --log-base 10.

Values are printed with 3 decimals by default (both per site and with -a).
Use --precision to print more (or fewer) decimals.

Characters are upper-cased before counting. The computation does not take into account
the following characters:
-> '*' and '.' (always)
-> '-' (if --remove-gaps is given)
-> every character that is not a standard nucleotide (A, C, G, T, U) or one of the
   20 standard amino acids, depending on the alignment alphabet: N, X, IUPAC
   ambiguity codes (R, Y, K, M, S, W, B, D, H, V), B, Z, J, etc.
   (if --ignore-ambiguities is given)

With --ignore-ambiguities, the alphabet matters:
-> If the alphabet is auto-detected as "both" (all characters are valid nucleotides and
   amino acids, e.g. a protein alignment made only of A, C, G, T, N...), the alignment is
   treated as nucleotides. Use --alphabet aa to force amino acids.
-> If the alphabet is unknown (characters that are neither nucleotides nor amino acids),
   an error is returned instead of a silently wrong result.

Without --ignore-ambiguities, N, X and ambiguity codes are counted as states of their
own, which may strongly inflate entropy on alignments with many undetermined positions.

If no character of a site is taken into account, then its entropy will be "NaN",
and it will not be taken into account in the average.
`,
	RunE: func(cmd *cobra.Command, args []string) (err error) {
		var aligns *align.AlignChannel
		var e float64
		var logfactor float64 = 1.0

		if entropyPrecision < 0 {
			err = fmt.Errorf("invalid --precision %d: must be >= 0", entropyPrecision)
			io.LogError(err)
			return
		}

		if entropyLogBase != "e" {
			var base float64
			if base, err = strconv.ParseFloat(entropyLogBase, 64); err != nil || base <= 0 || base == 1 {
				err = fmt.Errorf("invalid --log-base %q: must be \"e\" or a positive number different from 1", entropyLogBase)
				io.LogError(err)
				return
			}
			logfactor = math.Log(base)
		}

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
				if e, err = align.EntropyFiltered(i, entropyRemoveGaps, entropyRemoveAmbiguities); err != nil {
					io.LogError(err)
					return
				} else {
					e /= logfactor
					if entropyAverage {
						if !math.IsNaN(e) {
							avg += e
							total++
						}
					} else {
						fmt.Printf("%d\t%d\t%.*f\n", nb, i, entropyPrecision, e)
					}
				}
			}
			if entropyAverage {
				fmt.Printf("%d\t%.*f\n", nb, entropyPrecision, avg/float64(total))
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
	entropyCmd.PersistentFlags().BoolVar(&entropyRemoveAmbiguities, "ignore-ambiguities", false, "If true, then do not take into account N, X and other ambiguous/non standard characters in the computation")
	entropyCmd.PersistentFlags().StringVar(&entropyLogBase, "log-base", "e", "Base of the logarithm: \"e\" (nats, default), 2 (bits), 10, or any positive number")
	entropyCmd.PersistentFlags().IntVar(&entropyPrecision, "precision", 3, "Number of decimals of output entropy values (per site and average)")
}
