package cmd

import (
	"github.com/evolbioinfo/goalign/align"
	"github.com/evolbioinfo/goalign/io"
	"github.com/evolbioinfo/goalign/io/utils"
	"github.com/spf13/cobra"
)

var concatout string

// concatCmd represents the concat command
var concatCmd = &cobra.Command{
	Use:   "concat",
	Short: "Concatenates a set of alignments",
	Long: `Concatenates a set of alignments.
For example:

If format is phylip, it may contain several alignments in one file. 
Then we can concatenate all of them:
goalign concat -i align.phy

If format is Fasta, it is not possible, then you must give other alignments in the form:
goalign concat -i align.fasta others*.fasta

It is possible to give only otherfiles, without -i, by giving -i none
   goalign concat -i none align*.fasta
or goalign concat -i none -p align*.phy

`,
	RunE: func(cmd *cobra.Command, args []string) (err error) {
		var aligns *align.AlignChannel
		var alchan *align.AlignChannel
		var align align.Alignment = nil
		var f utils.StringWriterCloser

		if infile != "none" {
			if aligns, err = readalign(infile); err != nil {
				io.LogError(err)
				return
			}
			for al := range aligns.Achan {
				if align == nil {
					align = al
				} else {
					if err = align.Concat(al); err != nil {
						io.LogError(err)
						return
					}
				}
			}
			if aligns.Err != nil {
				err = aligns.Err
				io.LogError(err)
				return
			}
		}
		for _, otherfile := range args {
			if alchan, err = readalign(otherfile); err != nil {
				io.LogError(err)
				return
			}
			for al := range alchan.Achan {
				if align == nil {
					align = al
				} else {
					if err = align.Concat(al); err != nil {
						io.LogError(err)
						return
					}
				}
			}
			if alchan.Err != nil {
				err = alchan.Err
				io.LogError(err)
				return
			}
		}

		if f, err = utils.OpenWriteFile(concatout); err != nil {
			io.LogError(err)
			return
		}
		writeAlign(align, f)
		utils.CloseWriteFile(f, concatout)

		return
	},
}

func init() {
	RootCmd.AddCommand(concatCmd)
	concatCmd.PersistentFlags().StringVarP(&concatout, "output", "o", "stdout", "Alignment output file")
}
