package cmd

import (
	"bytes"
	"errors"
	"fmt"
	"os"

	"github.com/fredericlemoine/goalign/align"
	"github.com/fredericlemoine/goalign/io"
	"github.com/spf13/cobra"
)

var pssmlog bool
var pssmpseudocount float64
var pssmnorm int

// pssmCmd represents the pssm command
var pssmCmd = &cobra.Command{
	Use:   "pssm",
	Short: "Compute and print a Position specific scoring matrix",
	Long: `Compute and print a Position specific scoring matrix.

Different normalizations are possible with --normalization (-n):
	--normalization 0: None, means raw counts
	--normalization 1: By column frequency, i.e. frequency of nt/aa per site/column
	--normalization 2: By column frequency compared to alignment frequency: same as -n 1,
                           but divides by frequency of the nt/aa in the whole alignment
	--normalization 3: By column frequency compared to uniform frequency: same as -n 1, 
                           but divides by uniform frequency of the nt/aa (1/4 for nt, 1/20 for aa)
	--normalization 4: Normalization "Logo"

Possible to add pseudo counts (before normalization) with --pseudo-count (-c)

Possible to log2 transform the (normalized) value with --log (-l). Not taken into account with logo normalization
`,
	Run: func(cmd *cobra.Command, args []string) {
		switch pssmnorm {
		case align.PSSM_NORM_UNIF, align.PSSM_NORM_NONE, align.PSSM_NORM_FREQ, align.PSSM_NORM_DATA, align.PSSM_NORM_LOGO:
			for al := range rootaligns {
				if pssm, err := al.Pssm(pssmlog, pssmpseudocount, pssmnorm); err != nil {
					io.ExitWithMessage(err)
				} else {
					fmt.Fprintf(os.Stdout, printPSSM(al, pssm))
				}
			}

		default:
			io.ExitWithMessage(errors.New(fmt.Sprintf("Normlization does not exist: %d", pssmnorm)))
		}
	},
}

func init() {
	computeCmd.AddCommand(pssmCmd)
	pssmCmd.PersistentFlags().BoolVarP(&pssmlog, "log", "l", false, "(normalized) Values in log2")
	pssmCmd.PersistentFlags().Float64VarP(&pssmpseudocount, "pseudo-counts", "c", 0.0, "Value added to (normalized) counts")
	pssmCmd.PersistentFlags().IntVarP(&pssmnorm, "normalization", "n", 0, "Counts normalization")
}

func printPSSM(a align.Alignment, pssm map[rune][]float64) string {
	var buffer bytes.Buffer
	size := -1
	for _, c := range a.AlphabetCharacters() {
		v, ok := pssm[c]
		if !ok {
			io.ExitWithMessage(errors.New(fmt.Sprintf("Alphabet character %c is not in the pssm", c)))
		}
		buffer.WriteString(fmt.Sprintf("\t%c", c))

		if size == -1 {
			size = len(v)
		} else {
			if len(v) != size {
				io.ExitWithMessage(errors.New(fmt.Sprintf("Pssm has different sequence lengths for different characters")))
			}
		}
	}
	buffer.WriteString("\n")

	for i := 0; i < size; i++ {
		buffer.WriteString(fmt.Sprintf("%d", i+1))
		for _, c := range a.AlphabetCharacters() {
			v, _ := pssm[c]
			count := v[i]
			buffer.WriteString(fmt.Sprintf("\t%.3f", count))
		}
		buffer.WriteString("\n")
	}
	return buffer.String()
}
