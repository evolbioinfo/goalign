package cmd

import (
	"fmt"
	"os"

	"github.com/evolbioinfo/goalign/align"
	"github.com/evolbioinfo/goalign/io"
	"github.com/spf13/cobra"
)

var siteout string
var sitelength int
var sitenb int

// samplesitesCmd represents the samplesites command
var samplesitesCmd = &cobra.Command{
	Use:   "sites",
	Short: "Take a random subalignment",
	Long: `Take a random subalignment.

It take a random start position, and extract the alignment starting at that position
and with a given length.
`,
	RunE: func(cmd *cobra.Command, args []string) (err error) {
		var aligns *align.AlignChannel
		var f *os.File
		var subalign align.Alignment

		if aligns, err = readalign(infile); err != nil {
			io.LogError(err)
			return
		}

		al, _ := <-aligns.Achan
		if aligns.Err != nil {
			err = aligns.Err
			io.LogError(err)
			return
		}

		var name string = siteout
		var extension string = "fa"
		if rootphylip {
			extension = "phy"
		}
		for i := 0; i < sitenb; i++ {
			if sitenb > 1 {
				name = fmt.Sprintf("%s_%d.%s", siteout, i, extension)
			}
			if f, err = openWriteFile(name); err != nil {
				io.LogError(err)
				return
			}
			defer closeWriteFile(f, name)

			if subalign, err = al.RandSubAlign(sitelength); err != nil {
				io.LogError(err)
				return
			}
			writeAlign(subalign, f)
		}
		return
	},
}

func init() {
	sampleCmd.AddCommand(samplesitesCmd)
	samplesitesCmd.PersistentFlags().StringVarP(&siteout, "output", "o", "stdout", "Alignment output file")
	samplesitesCmd.PersistentFlags().IntVarP(&sitelength, "length", "l", 10, "Length of the random sub alignment")
	samplesitesCmd.PersistentFlags().IntVarP(&sitenb, "nsamples", "n", 1, "Number of samples to generate")
}
