package cmd

import (
	"os"

	"github.com/evolbioinfo/goalign/align"
	"github.com/evolbioinfo/goalign/io"
	"github.com/spf13/cobra"
)

var maskout string = "stdout"
var maskstart int
var masklength int
var maskunique bool
var maskrefseq string

// subseqCmd represents the subseq command
var maskCmd = &cobra.Command{
	Use:   "mask",
	Short: "Mask a sub alignment",
	Long: `Mask a part of the input alignment (replace by N|X)


It takes an alignment and replaces some characters by "N" or "X".

By default, it masks positions defined by start (0-based inclusive)
and length with N (nucleotide alignment) or X (amino-acid alignment).
If the length is after the end of the alignment, will stop at the 
end of the alignment.

For example:
goalign mask -p -i al.phy -s 9 -l 10

This will replace 10 positions with N|X from the 10th position.

If --ref-seq is specified, then coordinates are considered on the given reference sequence
without considering gaps. In that case, if the range of masked sites incorporates gaps in
the reference sequence, these gaps will also be masked in the output alignment.


If --unique is specified, 'goalign mask --unique' will replace characters that
are unique in their column (except GAPS) with N or X.
In this case, --length and --start are ignored.

The output format is the same than input format.
`,
	RunE: func(cmd *cobra.Command, args []string) (err error) {
		var aligns *align.AlignChannel
		var f *os.File

		if aligns, err = readalign(infile); err != nil {
			io.LogError(err)
			return
		}
		if f, err = openWriteFile(maskout); err != nil {
			io.LogError(err)
			return
		}

		refseq := cmd.Flags().Changed("ref-seq")

		for al := range aligns.Achan {
			if aligns.Err != nil {
				err = aligns.Err
				io.LogError(err)
			}
			if maskunique {
				if err = al.MaskUnique(); err != nil {
					io.LogError(err)
					return
				}
			} else {
				start := maskstart
				length := masklength
				if refseq {
					start, length, err = al.RefCoordinates(maskrefseq, start, length)
				}
				if err = al.Mask(start, length); err != nil {
					io.LogError(err)
					return
				}
			}
			writeAlign(al, f)
		}
		f.Close()

		return
	},
}

func init() {
	RootCmd.AddCommand(maskCmd)
	maskCmd.PersistentFlags().StringVarP(&maskout, "output", "o", "stdout", "Alignment output file")
	maskCmd.PersistentFlags().IntVarP(&maskstart, "start", "s", 0, "Start position (0-based inclusive)")
	maskCmd.PersistentFlags().IntVarP(&masklength, "length", "l", 10, "Length of the sub alignment")
	maskCmd.PersistentFlags().StringVar(&maskrefseq, "ref-seq", "none", "Coordinates are considered wrt. to the given reference sequence (no effect with --unique)")
	maskCmd.PersistentFlags().BoolVar(&maskunique, "unique", false, "If given, then masks characters that are unique in their columns (start and length are ignored)")
}
