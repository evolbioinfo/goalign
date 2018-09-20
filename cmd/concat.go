package cmd

import (
	"github.com/spf13/cobra"

	"github.com/fredericlemoine/goalign/align"
	"github.com/fredericlemoine/goalign/io"
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
	Run: func(cmd *cobra.Command, args []string) {
		var align align.Alignment = nil
		var err error
		if infile != "none" {
			aligns := readalign(infile)
			for al := range aligns.Achan {
				if align == nil {
					align = al
				} else {
					err = align.Concat(al)
					if err != nil {
						io.ExitWithMessage(err)
					}
				}
			}
			if aligns.Err != nil {
				io.ExitWithMessage(aligns.Err)
			}
		}
		for _, otherfile := range args {
			alchan := readalign(otherfile)
			for al := range alchan.Achan {
				if align == nil {
					align = al
				} else {
					err = align.Concat(al)
					if err != nil {
						io.ExitWithMessage(err)
					}
				}
			}
			if alchan.Err != nil {
				io.ExitWithMessage(alchan.Err)
			}
		}

		f := openWriteFile(concatout)
		writeAlign(align, f)
		f.Close()
	},
}

func init() {
	RootCmd.AddCommand(concatCmd)
	concatCmd.PersistentFlags().StringVarP(&concatout, "output", "o", "stdout", "Alignment output file")

	// Here you will define your flags and configuration settings.

	// Cobra supports Persistent Flags which will work for this command
	// and all subcommands, e.g.:
	// concatCmd.PersistentFlags().String("foo", "", "A help for foo")

	// Cobra supports local flags which will only run when this command
	// is called directly, e.g.:
	// concatCmd.Flags().BoolP("toggle", "t", false, "Help message for toggle")

}
