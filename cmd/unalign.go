package cmd

import (
	"fmt"

	"github.com/evolbioinfo/goalign/align"
	"github.com/evolbioinfo/goalign/io"
	"github.com/evolbioinfo/goalign/io/utils"
	"github.com/spf13/cobra"
)

var unalignOutput string

// unalignCmd represents the unalign command
var unalignCmd = &cobra.Command{
	Use:   "unalign",
	Short: "Unaligns input alignment",
	Long: `Unaligns an input alignment, by removing indels.

The output is in Fasta format, whatever the input format is (Fasta or Phylip).

Output sequences are free from all indel charachers "-".

As there may be several alignments in the input alignment (phylip format),
output files are prefixed with argument given to "--output-prefix", and 
suffixed with an index and the extension ".fa".

If --output-prefix is set to "-" or "stdout", all sequences are printed on stdout

Example:

goalign unalign -i align.ph -p -o seq_

If align contains 3 alignments, this will generate 3 files:
* seq_000001.fa
* seq_000002.fa
* seq_000003.fa
`,
	RunE: func(cmd *cobra.Command, args []string) (err error) {
		var aligns *align.AlignChannel
		var f utils.StringWriterCloser

		if aligns, err = readalign(infile); err != nil {
			io.LogError(err)
			return
		}

		i := 1
		filename := unalignOutput
		for al := range aligns.Achan {
			if filename != "stdout" && filename != "-" {
				filename = fmt.Sprintf("%s_%.6d.fa", unalignOutput, i)
			}
			if f, err = utils.OpenWriteFile(filename); err != nil {
				io.LogError(err)
				return
			}
			writeSequences(al.Unalign(), f)
			utils.CloseWriteFile(f, filename)
			i++
		}

		if aligns.Err != nil {
			err = aligns.Err
			io.LogError(err)
		}
		return
	},
}

func init() {
	RootCmd.AddCommand(unalignCmd)
	unalignCmd.PersistentFlags().StringVarP(&unalignOutput, "output-prefix", "o", "stdout", "Unaligned alignment output file prefix")
}
