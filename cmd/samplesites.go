package cmd

import (
	"fmt"

	"github.com/fredericlemoine/goalign/io"
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
	Run: func(cmd *cobra.Command, args []string) {
		aligns := readalign(infile)
		al, _ := <-aligns.Achan
		if aligns.Err != nil {
			io.ExitWithMessage(aligns.Err)
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
			out := openWriteFile(name)
			subalign, err := al.RandSubAlign(sitelength)
			if err != nil {
				io.ExitWithMessage(err)
			}
			writeAlign(subalign, out)
			out.Close()
		}
	},
}

func init() {
	sampleCmd.AddCommand(samplesitesCmd)
	samplesitesCmd.PersistentFlags().StringVarP(&siteout, "output", "o", "stdout", "Alignment output file")
	samplesitesCmd.PersistentFlags().IntVarP(&sitelength, "length", "l", 10, "Length of the random sub alignment")
	samplesitesCmd.PersistentFlags().IntVarP(&sitenb, "nsamples", "n", 1, "Number of samples to generate")
}
