package cmd

import (
	"bytes"
	"fmt"
	"os"

	"github.com/evolbioinfo/goalign/align"
	"github.com/evolbioinfo/goalign/io"
	"github.com/spf13/cobra"
)

var pssmlog bool
var pssmpseudocount float64
var pssmnorm int

// pssmCmd represents the pssm command
var pssmCmd = &cobra.Command{
	Use:   "pssm",
	Short: "Computes and prints a Position specific scoring matrix",
	Long: `Computes and prints a Position specific scoring matrix.

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
	RunE: func(cmd *cobra.Command, args []string) (err error) {
		var aligns *align.AlignChannel
		var pssm map[uint8][]float64
		var pssmstring string

		if aligns, err = readalign(infile); err != nil {
			io.LogError(err)
			return
		}

		switch pssmnorm {
		case align.PSSM_NORM_UNIF, align.PSSM_NORM_NONE, align.PSSM_NORM_FREQ, align.PSSM_NORM_DATA, align.PSSM_NORM_LOGO:
			for al := range aligns.Achan {
				if pssm, err = al.Pssm(pssmlog, pssmpseudocount, pssmnorm); err != nil {
					io.LogError(err)
					return
				}
				if pssmstring, err = printPSSM(al, pssm); err != nil {
					io.LogError(err)
					return
				}
				fmt.Fprintf(os.Stdout, pssmstring)
			}

		default:
			err = fmt.Errorf("Normlization does not exist: %d", pssmnorm)
			io.LogError(err)
			return
		}

		if aligns.Err != nil {
			err = aligns.Err
			io.LogError(err)
		}

		return
	},
}

func init() {
	computeCmd.AddCommand(pssmCmd)
	pssmCmd.PersistentFlags().BoolVarP(&pssmlog, "log", "l", false, "(normalized) Values in log2")
	pssmCmd.PersistentFlags().Float64VarP(&pssmpseudocount, "pseudo-counts", "c", 0.0, "Value added to (normalized) counts")
	pssmCmd.PersistentFlags().IntVarP(&pssmnorm, "normalization", "n", 0, "Counts normalization")
}

func printPSSM(a align.Alignment, pssm map[uint8][]float64) (pssmstring string, err error) {
	var buffer bytes.Buffer
	size := -1
	for _, c := range a.AlphabetCharacters() {
		v, ok := pssm[c]
		if !ok {
			err = fmt.Errorf("Alphabet character %c is not in the pssm", c)
			return
		}
		buffer.WriteString(fmt.Sprintf("\t%c", c))

		if size == -1 {
			size = len(v)
		} else {
			if len(v) != size {
				err = fmt.Errorf("Pssm has different sequence lengths for different characters")
				return
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
	pssmstring = buffer.String()
	return
}
