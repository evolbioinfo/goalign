package cmd

import (
	"github.com/evolbioinfo/goalign/align"
	"github.com/evolbioinfo/goalign/io"
	"github.com/evolbioinfo/goalign/io/utils"
	"github.com/spf13/cobra"
)

var maskout string = "stdout"
var maskstart int
var masklength int
var maskatmost int
var maskunique bool
var maskrefseq string
var maskreplace string

// subseqCmd represents the subseq command
var maskCmd = &cobra.Command{
	Use:   "mask",
	Short: "Mask a sub alignment",
	Long: `Mask a part of the input alignment (replace by N|X by default)


It takes an alignment and replaces some characters by "N" or "X" by default.

By default, it masks positions defined by start (0-based inclusive)
and length with N (nucleotide alignment) or X (amino-acid alignment) by default.
If the length is after the end of the alignment, will stop at the 
end of the alignment.

For example:
goalign mask -p -i al.phy -s 9 -l 10

This will replace 10 positions with N|X from the 10th position.

If --ref-seq is specified, then coordinates are considered on the given reference sequence
without considering gaps. In that case, if the range of masked sites incorporates gaps in
the reference sequence, these gaps will also be masked in the output alignment.

If --unique is specified, 'goalign mask --unique' will replace characters that
are unique (defined by --at-most option) in their column (except GAPS) with N or X. 
In that case, if --ref-seq option is given, then a unique character is masked if:
    1) It is different from the given reference sequence
    2) Or the reference is a GAP
In this case, --length and --start are ignored.

Option --replace defines the replacement character. If --replace is "" (default), then, 
masked characters are replaced by "N" or "X" depending on the alphabet. 
Orherwise:
  1) if --replace is AMBIG: just like ""
  2) if --replace is MAJ: Masked characters are replaced by the most frequent character of the column 
     (without considering the reference sequence when --ref-seq is given)
  3) if --replace is GAP: Replacing character is a GAP

The output format is the same than input format.
`,
	RunE: func(cmd *cobra.Command, args []string) (err error) {
		var aligns *align.AlignChannel
		var f utils.StringWriterCloser

		if aligns, err = readalign(infile); err != nil {
			io.LogError(err)
			return
		}
		if f, err = utils.OpenWriteFile(maskout); err != nil {
			io.LogError(err)
			return
		}

		refseq := cmd.Flags().Changed("ref-seq")

		if !refseq {
			maskrefseq = ""
		}

		for al := range aligns.Achan {
			if aligns.Err != nil {
				err = aligns.Err
				io.LogError(err)
			}
			if maskunique {
				if err = al.MaskOccurences(maskrefseq, maskatmost, maskreplace); err != nil {
					io.LogError(err)
					return
				}
			} else {
				start := maskstart
				length := masklength
				if refseq {
					if start, length, err = al.RefCoordinates(maskrefseq, start, length); err != nil {
						io.LogError(err)
						return
					}
				}
				if err = al.Mask(start, length, maskreplace); err != nil {
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
	maskCmd.PersistentFlags().StringVar(&maskrefseq, "ref-seq", "none", "Coordinates are considered wrt. to the given reference sequence (with --unique, it masks unique characters that are different from the reference sequence)")
	maskCmd.PersistentFlags().BoolVar(&maskunique, "unique", false, "If given, then masks characters that are unique (defined with --at-most) in their columns (start and length are ignored)")
	maskCmd.PersistentFlags().StringVar(&maskreplace, "replace", "AMBIG", "Replacement character. If AMBIG: N or X (depending on alphabet), if GAP: -, if MAJ: the main character of the column, or can be any other character")
	maskCmd.PersistentFlags().IntVar(&maskatmost, "at-most", 1, "The number of occurences that defines the uniqueness of the characher in the column (only used with --unique)")
}
