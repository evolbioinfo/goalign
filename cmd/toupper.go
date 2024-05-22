package cmd

import (
	"github.com/evolbioinfo/goalign/align"
	"github.com/evolbioinfo/goalign/io"
	"github.com/evolbioinfo/goalign/io/utils"
	"github.com/spf13/cobra"
)

var toupperOutput string

var toupperCmd = &cobra.Command{
	Use:   "toupper",
	Short: "Replace lower case characters by upper case characters",
	Long: `Replace lower case characters by upper case characters.
`,
	RunE: func(cmd *cobra.Command, args []string) (err error) {
		var f utils.StringWriterCloser

		if f, err = utils.OpenWriteFile(toupperOutput); err != nil {
			io.LogError(err)
			return
		}
		defer utils.CloseWriteFile(f, toupperOutput)

		if unaligned {
			var seqs align.SeqBag

			if seqs, err = readsequences(infile); err != nil {
				io.LogError(err)
				return
			}
			seqs.ToUpper()
			writeSequences(seqs, f)
		} else {
			var aligns *align.AlignChannel
			var al align.Alignment

			if aligns, err = readalign(infile); err != nil {
				io.LogError(err)
				return
			}
			for al = range aligns.Achan {
				al.ToUpper()
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
	RootCmd.AddCommand(toupperCmd)
	toupperCmd.PersistentFlags().StringVarP(&toupperOutput, "output", "o", "stdout", "Output upper case alignment file")
	toupperCmd.PersistentFlags().BoolVar(&unaligned, "unaligned", false, "Considers sequences as unaligned and format fasta (phylip, nexus,... options are ignored)")
}
