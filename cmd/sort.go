package cmd

import (
	"os"

	"github.com/evolbioinfo/goalign/align"
	"github.com/evolbioinfo/goalign/io"
	"github.com/spf13/cobra"
)

var sortOutput string

// reformatCmd represents the reformat command
var sortCmd = &cobra.Command{
	Use:   "sort",
	Short: "sorts input alignment by sequence name",
	Long: `sorts input algignment by sequence name.
`,
	RunE: func(cmd *cobra.Command, args []string) (err error) {
		var f *os.File

		if f, err = openWriteFile(sortOutput); err != nil {
			io.LogError(err)
			return
		}
		defer closeWriteFile(f, sortOutput)

		if unaligned {
			var seqs align.SeqBag
			if seqs, err = readsequences(infile); err != nil {
				io.LogError(err)
				return
			}
			seqs.Sort()
			writeSequences(seqs, f)
		} else {
			var aligns *align.AlignChannel

			if aligns, err = readalign(infile); err != nil {
				io.LogError(err)
				return
			}

			for al := range aligns.Achan {
				al.Sort()
				writeAlign(al, f)
			}

			if aligns.Err != nil {
				err = aligns.Err
				io.LogError(err)
			}
		}
		return
	},
}

func init() {
	sortCmd.PersistentFlags().BoolVar(&unaligned, "unaligned", false, "Considers sequences as unaligned and format fasta (phylip, nexus,... options are ignored)")
	sortCmd.PersistentFlags().StringVarP(&sortOutput, "output", "o", "stdout", "Sorted alignment output file")
	RootCmd.AddCommand(sortCmd)
}
