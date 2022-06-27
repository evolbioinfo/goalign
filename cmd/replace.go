package cmd

import (
	"errors"

	"github.com/evolbioinfo/goalign/align"
	"github.com/evolbioinfo/goalign/io"
	"github.com/evolbioinfo/goalign/io/utils"
	"github.com/spf13/cobra"
)

var replaceOutput string
var replaceRegexp bool
var replaceOld, replaceNew string

// renameCmd represents the rename command
var replaceCmd = &cobra.Command{
	Use:   "replace",
	Short: "Replace characters in sequences of the input alignment (possible with a regex)",
	Long: `Replace characters in sequences of the input alignment (possible with a regex).
If the replacement changes sequence length, then returns an error.
`,
	RunE: func(cmd *cobra.Command, args []string) (err error) {
		var aligns *align.AlignChannel
		var f utils.StringWriterCloser
		var seqs align.SeqBag

		if !cmd.Flags().Changed("old") || !cmd.Flags().Changed("new") {
			err = errors.New("--old and --new must be specified")
			return
		}

		if f, err = utils.OpenWriteFile(replaceOutput); err != nil {
			io.LogError(err)
			return
		}
		defer utils.CloseWriteFile(f, replaceOutput)

		if unaligned {
			if seqs, err = readsequences(infile); err != nil {
				io.LogError(err)
				return
			}
			if err = seqs.Replace(replaceOld, replaceNew, replaceRegexp); err != nil {
				io.LogError(err)
				return
			}
			writeSequences(seqs, f)
		} else {
			if aligns, err = readalign(infile); err != nil {
				io.LogError(err)
				return
			}
			for al := range aligns.Achan {
				if err = al.Replace(replaceOld, replaceNew, replaceRegexp); err != nil {
					io.LogError(err)
					return
				}
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
	RootCmd.AddCommand(replaceCmd)

	replaceCmd.PersistentFlags().StringVarP(&replaceOutput, "output", "o", "stdout", "Output alignment file")
	replaceCmd.PersistentFlags().BoolVarP(&replaceRegexp, "regexp", "e", false, "Considers Replace alignment using regexp")
	replaceCmd.PersistentFlags().StringVarP(&replaceOld, "old", "s", "none", "String to replace in the sequences")
	replaceCmd.PersistentFlags().StringVarP(&replaceNew, "new", "n", "none", "New string that will replace old string in sequences")
	replaceCmd.PersistentFlags().BoolVar(&unaligned, "unaligned", false, "Considers input sequences as unaligned and fasta format (phylip, nexus,... options are ignored)")
}
