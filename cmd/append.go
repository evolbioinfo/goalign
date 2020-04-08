package cmd

import (
	"os"

	"github.com/evolbioinfo/goalign/align"
	"github.com/evolbioinfo/goalign/io"
	"github.com/spf13/cobra"
)

var appendout string

// appendCmd represents the append command
var appendCmd = &cobra.Command{
	Use:   "append",
	Short: "Append alignments to an input alignment",
	Long: `Append alignments to an input alignment.

This commands adds the sequences of a set of alignments to a reference alignement
specified by -i.

If sequences do not have the same length than the reference alignment, then returns an error.

If format is phylip, it may contain several alignments in one file. 
Then we can append all of them at once:
goalign append -i refalign.phy aligns.phy

If format is Fasta, several alignments may be given in the form:
goalign append -i align.fasta others*.fasta

`,
	RunE: func(cmd *cobra.Command, args []string) (err error) {
		var compAligns *align.AlignChannel
		var refAligns *align.AlignChannel = nil
		var refAlign align.Alignment = nil

		var f *os.File

		if infile != "none" {
			if refAligns, err = readalign(infile); err != nil {
				io.LogError(err)
				return
			}
			for al := range refAligns.Achan {
				if refAlign == nil {
					refAlign = al
				} else {
					if err = refAlign.Append(al); err != nil {
						io.LogError(err)
						return
					}
				}
				if refAligns.Err != nil {
					err = refAligns.Err
					io.LogError(err)
					return
				}
			}
		}

		for _, otherfile := range args {
			if compAligns, err = readalign(otherfile); err != nil {
				io.LogError(err)
				return
			}
			for al := range compAligns.Achan {
				if err = refAlign.Append(al); err != nil {
					io.LogError(err)
					return
				}
			}
			if compAligns.Err != nil {
				err = compAligns.Err
				io.LogError(err)
				return
			}
		}

		if f, err = openWriteFile(appendout); err != nil {
			io.LogError(err)
			return
		}
		writeAlign(refAlign, f)
		closeWriteFile(f, appendout)

		return
	},
}

func init() {
	RootCmd.AddCommand(appendCmd)
	appendCmd.PersistentFlags().StringVarP(&appendout, "output", "o", "stdout", "Alignment output file")
}
