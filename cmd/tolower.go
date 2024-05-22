package cmd

import (
	"github.com/evolbioinfo/goalign/align"
	"github.com/evolbioinfo/goalign/io"
	"github.com/evolbioinfo/goalign/io/utils"
	"github.com/spf13/cobra"
)

var tolowerOutput string

var tolowerCmd = &cobra.Command{
	Use:   "tolower",
	Short: "Replace upper case characters by lower case characters",
	Long: `Replace upper case characters by lower case characters.
`,
	RunE: func(cmd *cobra.Command, args []string) (err error) {
		var f utils.StringWriterCloser

		if f, err = utils.OpenWriteFile(tolowerOutput); err != nil {
			io.LogError(err)
			return
		}
		defer utils.CloseWriteFile(f, tolowerOutput)

		if unaligned {
			var seqs align.SeqBag

			if seqs, err = readsequences(infile); err != nil {
				io.LogError(err)
				return
			}
			seqs.ToLower()
			writeSequences(seqs, f)
		} else {
			var aligns *align.AlignChannel
			var al align.Alignment

			if aligns, err = readalign(infile); err != nil {
				io.LogError(err)
				return
			}
			for al = range aligns.Achan {
				al.ToLower()
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
	RootCmd.AddCommand(tolowerCmd)
	tolowerCmd.PersistentFlags().StringVarP(&tolowerOutput, "output", "o", "stdout", "Output lower case alignment file")
	tolowerCmd.PersistentFlags().BoolVar(&unaligned, "unaligned", false, "Considers sequences as unaligned and format fasta (phylip, nexus,... options are ignored)")
}
