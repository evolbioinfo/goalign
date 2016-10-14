package cmd

import (
	"fmt"

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
	Run: func(cmd *cobra.Command, args []string) {
		i := 1
		filename := unalignOutput
		for al := range rootaligns {
			if filename != "stdout" && filename != "-" {
				filename = fmt.Sprintf("%s_%.6d.fa", unalignOutput, i)
			}
			f := openWriteFile(filename)
			writeUnAlignFasta(al, f)
			if filename != "stdout" && filename != "-" {
				f.Close()
			}
		}
	},
}

func init() {
	RootCmd.AddCommand(unalignCmd)
	unalignCmd.PersistentFlags().StringVarP(&unalignOutput, "output-prefix", "o", "stdout", "Unaligned alignment output file prefix")
}
