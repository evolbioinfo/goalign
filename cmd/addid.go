package cmd

import (
	"github.com/evolbioinfo/goalign/align"
	"github.com/evolbioinfo/goalign/io"
	"github.com/evolbioinfo/goalign/io/utils"
	"github.com/spf13/cobra"
)

var addIdOutput string
var addIdName string
var addIdRight bool

// addidCmd represents the addid command
var addidCmd = &cobra.Command{
	Use:   "addid",
	Short: "Adds a string to each sequence identifier of the input alignment",
	Long: `This command adds an indentifier (string) to all sequences of an input alignment. 

The string may be added to the left or to the right of each sequence name.
`,
	RunE: func(cmd *cobra.Command, args []string) (err error) {
		var f utils.StringWriterCloser

		if f, err = utils.OpenWriteFile(addIdOutput); err != nil {
			io.LogError(err)
			return
		}
		defer utils.CloseWriteFile(f, addIdOutput)

		if unaligned {
			var seqs align.SeqBag

			if seqs, err = readsequences(infile); err != nil {
				io.LogError(err)
				return
			}
			seqs.AppendSeqIdentifier(addIdName, addIdRight)
			writeSequences(seqs, f)
		} else {

			var aligns *align.AlignChannel
			if aligns, err = readalign(infile); err != nil {
				io.LogError(err)
				return
			}
			for al := range aligns.Achan {
				al.AppendSeqIdentifier(addIdName, addIdRight)
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
	RootCmd.AddCommand(addidCmd)
	addidCmd.PersistentFlags().StringVarP(&addIdOutput, "out-align", "o", "stdout", "Renamed alignment output file")
	addidCmd.PersistentFlags().StringVarP(&addIdName, "name", "n", "none", "String to add to sequence names")
	addidCmd.PersistentFlags().BoolVarP(&addIdRight, "right", "r", false, "Adds the String on the right of sequence names (otherwise, adds to left)")
	addidCmd.PersistentFlags().BoolVar(&unaligned, "unaligned", false, "Considers sequences as unaligned and format fasta (phylip, nexus,... options are ignored)")
}
