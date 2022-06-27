package cmd

import (
	"fmt"

	"github.com/evolbioinfo/goalign/align"
	"github.com/evolbioinfo/goalign/io"
	"github.com/evolbioinfo/goalign/io/utils"
	"github.com/spf13/cobra"
)

var compressOutput string
var compressWeightOutput string

var compressCmd = &cobra.Command{
	Use:   "compress",
	Short: "Removes identical patterns/sites from an input alignment",
	Long: `Removes identical patterns/sites from an input alignment

And prints in the weight file the number of occurence of each pattern

Example: 

ali.phy
1 GGGGGGGGGGGGGGGGGGGG
2 TTTTTTTTTTTTTTTTTTTT
3 GGGGGGGGGGCCCCCCCCCC
4 AAAAAAAAAAAAAAAAAAAA

goalign compress -i  ali.phy will produce:
1 GG
2 TT
3 GC
4 AA

and weight file:
10
10
`,
	RunE: func(cmd *cobra.Command, args []string) (err error) {
		var aligns *align.AlignChannel
		var f, wf utils.StringWriterCloser

		if aligns, err = readalign(infile); err != nil {
			io.LogError(err)
			return
		}
		if f, err = utils.OpenWriteFile(compressOutput); err != nil {
			io.LogError(err)
			return
		}
		defer utils.CloseWriteFile(f, compressOutput)

		if wf, err = utils.OpenWriteFile(compressWeightOutput); err != nil {
			io.LogError(err)
			return
		}
		defer utils.CloseWriteFile(wf, compressWeightOutput)

		for al := range aligns.Achan {
			var w []int
			if w = al.Compress(); err != nil {
				io.LogError(err)
				return
			} else {
				writeAlign(al, f)
				writeWeights(w, wf)
			}
		}

		if aligns.Err != nil {
			err = aligns.Err
			io.LogError(err)
		}
		return
	},
}

func init() {
	compressCmd.PersistentFlags().StringVarP(&compressOutput, "output", "o", "stdout", "Compressed output alignment file")
	compressCmd.PersistentFlags().StringVar(&compressWeightOutput, "weight-out", "none", "Pattern weight output file")
	RootCmd.AddCommand(compressCmd)
}

func writeWeights(weights []int, f utils.StringWriterCloser) {
	for _, w := range weights {
		fmt.Fprintf(f, "%d\n", w)
	}
}
