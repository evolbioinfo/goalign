package cmd

import (
	"github.com/evolbioinfo/goalign/align"
	"github.com/evolbioinfo/goalign/io"
	"github.com/evolbioinfo/goalign/io/utils"
	"github.com/spf13/cobra"
)

var orfOutput string
var orfreverse bool

// translateCmd represents the addid command
var orfCmd = &cobra.Command{
	Use:   "orf",
	Short: "Find the longest orf in all given sequences in forward strand",
	Long: `Find the longest orf in all given sequences in forward strand.

If input sequences are not nucleotidic, then returns an error.
If input sequences are aligned (contain '-'), then they are unaligned first.

Output is in fasta format.
`,
	RunE: func(cmd *cobra.Command, args []string) (err error) {
		var f utils.StringWriterCloser
		var reforf align.SeqBag
		var inseqs align.SeqBag
		var orf align.Sequence

		if f, err = utils.OpenWriteFile(orfOutput); err != nil {
			io.LogError(err)
			return
		}
		defer utils.CloseWriteFile(f, orfOutput)

		if inseqs, err = readsequences(infile); err != nil {
			io.LogError(err)
			return
		}

		inseqs = inseqs.Unalign()

		if orf, err = inseqs.LongestORF(orfreverse); err != nil {
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
	orfCmd.PersistentFlags().BoolVar(&orfreverse, "reverse", false, "Search for the longest ORF ALSO in the reverse strand")
}
