package cmd

import (
	"os"

	"github.com/fredericlemoine/goalign/align"
	"github.com/fredericlemoine/goalign/io"
	"github.com/spf13/cobra"
)

var orfOutput string

// translateCmd represents the addid command
var orfCmd = &cobra.Command{
	Use:   "orf",
	Short: "Find the longest orf in all given sequences in forward strand",
	Long: `Find the longest orf in all given sequences in forward strand.

If input sequences are not nucleotidic, then returns an error.
If input sequences are aligned (contain '-'), then they are unaligned first.
`,
	RunE: func(cmd *cobra.Command, args []string) (err error) {
		var f *os.File
		var reforf align.SeqBag
		var inseqs align.SeqBag
		var orf align.Sequence

		if f, err = openWriteFile(orfOutput); err != nil {
			io.LogError(err)
			return
		}
		defer closeWriteFile(f, orfOutput)

		if inseqs, err = readsequences(infile); err != nil {
			io.LogError(err)
			return
		}

		inseqs = inseqs.Unalign()

		if orf, err = inseqs.LongestORF(); err != nil {
			io.LogError(err)
			return
		}
		reforf = align.NewSeqBag(align.UNKNOWN)
		reforf.AddSequenceChar(orf.Name(), orf.SequenceChar(), orf.Comment())
		reforf.AutoAlphabet()
		writeSequences(reforf, f)

		return
	},
}

func init() {
	RootCmd.AddCommand(orfCmd)
	orfCmd.PersistentFlags().StringVarP(&orfOutput, "output", "o", "stdout", "ORF Output Fasta File")
}
