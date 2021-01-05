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
var siteconsecutive bool

// samplesitesCmd represents the samplesites command
var samplesitesCmd = &cobra.Command{
	Use:   "sites",
	Short: "Take a random subalignment",
	Long: `Take a random subalignment.

	goalign sample sites extracts a subalignment of given length from the input alignment.
	* If --consecutive is true, then a start position is randomly chosen, and the next "length" 
	positions are extracted.
	* Otherwise, if consecutive is false, then "length" positions are sampled without replacement 
	from the original alignment (any order).

	If more than 2 samples are requested (-n x, x>1), then alignemnts are written in output files
	with suffix _0.ext, _1.ext, etc. (with ext in [ph,fa,clustal,nx] ).The only exception is when
	format is phylip and output file name is stdout	or -, then output alignemnts are written on stdout.
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
		var extension string = alignExtension()

		for i := 0; i < sitenb; i++ {
			if sitenb > 1 && !(rootphylip && (siteout == "stdout" || siteout == "-")) {
				name = fmt.Sprintf("%s_%d.%s", siteout, i, extension)
			}
			if f, err = openWriteFile(name); err != nil {
				io.LogError(err)
				return
			}
			defer closeWriteFile(f, name)

			if subalign, err = al.RandSubAlign(sitelength, siteconsecutive); err != nil {
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
	samplesitesCmd.PersistentFlags().BoolVar(&siteconsecutive, "consecutive", true, "If sampled sites are consecutive (inactivate with --consecutive=false)")
}
