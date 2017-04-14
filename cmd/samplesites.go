package cmd

import (
	"fmt"
	"math/rand"
	"time"

	"github.com/fredericlemoine/goalign/io"
	"github.com/spf13/cobra"
)

var siteseed int64
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
		rand.Seed(siteseed)

		al := <-rootaligns
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
	samplesitesCmd.PersistentFlags().Int64VarP(&siteseed, "seed", "s", time.Now().UTC().UnixNano(), "Initial Random Seed")
	samplesitesCmd.PersistentFlags().StringVarP(&siteout, "output", "o", "stdout", "Alignment output file")
	samplesitesCmd.PersistentFlags().IntVarP(&sitelength, "length", "l", 10, "Length of the random sub alignment")
	samplesitesCmd.PersistentFlags().IntVarP(&sitenb, "nsamples", "n", 1, "Number of samples to generate")
}
