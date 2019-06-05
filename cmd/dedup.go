package cmd

import (
	"os"

	"github.com/evolbioinfo/goalign/align"
	"github.com/evolbioinfo/goalign/io"
	"github.com/spf13/cobra"
)

var dedupOutput string
var dedupLogOutput string

// dedupCmd represents the dedup command
var dedupCmd = &cobra.Command{
	Use:   "dedup",
	Short: "Deduplicate sequences that have the same sequence",
	Long: `Deduplicate sequences that have the same sequence

The name of the first sequence is kept

Example: 

ali.phy
1 AAAAAA
2 CCCCCC
3 GGGGGG
4 GGGGGG

goalign dedup -i ali.phy will produce:

1 AAAAAA
2 CCCCCC
3 GGGGGG

if -l is specified, then identical sequences are printed in the given file
with the following format:

seq1,seq2
seq3,seq4

This means that seq1 is identical to seq2 and seq3 is identical to seq4.

`,
	RunE: func(cmd *cobra.Command, args []string) (err error) {
		var aligns *align.AlignChannel
		var f, l *os.File
		var id [][]string

		if aligns, err = readalign(infile); err != nil {
			io.LogError(err)
			return
		}
		if f, err = openWriteFile(dedupOutput); err != nil {
			io.LogError(err)
			return
		}
		defer closeWriteFile(f, dedupOutput)

		if l, err = openWriteFile(dedupLogOutput); err != nil {
			io.LogError(err)
			return
		}
		defer closeWriteFile(l, dedupLogOutput)

		for al := range aligns.Achan {
			if id, err = al.Deduplicate(); err != nil {
				io.LogError(err)
				return
			} else {
				writeAlign(al, f)
				writeIdentical(id, l)
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
	RootCmd.AddCommand(dedupCmd)
	dedupCmd.PersistentFlags().StringVarP(&dedupOutput, "output", "o", "stdout", "Deduplicated output alignment file")
	dedupCmd.PersistentFlags().StringVarP(&dedupLogOutput, "log", "l", "none", "Deduplicated output log file")
}

func writeIdentical(id [][]string, logfile *os.File) {
	for _, s := range id {
		for i, name := range s {
			if i > 0 {
				logfile.WriteString(",")
			}
			logfile.WriteString(name)
		}
		logfile.WriteString("\n")
	}
}
