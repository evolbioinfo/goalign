package cmd

import (
	"github.com/evolbioinfo/goalign/align"
	"github.com/evolbioinfo/goalign/io"
	"github.com/evolbioinfo/goalign/io/utils"
	"github.com/spf13/cobra"
)

var revCompOutput string

// revCompCmd
var revCompCmd = &cobra.Command{
	Use:   "revcomp",
	Short: "Prints the reverse complement of all sequences",
	Long: `Prints the reverse complement of all sequences of the alignement.
	
	Sequences may be unaligned (--unaligned option).

If the input alignment is not nucleotides, then returns an error.

IUPAC codes are taken into account for the reverse complement.

If sequence names are given in the command line (e.g. goalign revcomp -i al.fasta s1 s2 s3), 
only given sequences are reverse-complemented, if they exist in the alignment.
`,
	RunE: func(cmd *cobra.Command, args []string) (err error) {
		var f utils.StringWriterCloser

		if f, err = utils.OpenWriteFile(revCompOutput); err != nil {
			io.LogError(err)
			return
		}
		defer utils.CloseWriteFile(f, revCompOutput)

		if unaligned {
			var seqs align.SeqBag

			if seqs, err = readsequences(infile); err != nil {
				io.LogError(err)
				return
			}

			// If names are given as arguments, we reverse complement only these sequences
			if len(args) > 0 {
				if err = seqs.ReverseComplementSequences(args...); err != nil {
					io.LogError(err)
					return
				}
			} else {
				if err = seqs.ReverseComplement(); err != nil {
					io.LogError(err)
					return
				}
			}
			writeSequences(seqs, f)
		} else {
			var aligns *align.AlignChannel
			var al align.Alignment

			if aligns, err = readalign(infile); err != nil {
				io.LogError(err)
				return
			}
			for al = range aligns.Achan {
				// If names are given as arguments, we reverse complement only these sequences
				if len(args) > 0 {
					if err = al.ReverseComplementSequences(args...); err != nil {
						io.LogError(err)
						return
					}
				} else {
					if err = al.ReverseComplement(); err != nil {
						io.LogError(err)
						return
					}
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
	RootCmd.AddCommand(revCompCmd)
	revCompCmd.PersistentFlags().StringVarP(&revCompOutput, "output", "o", "stdout", "Output reverse complement alignment file")
	revCompCmd.PersistentFlags().BoolVar(&unaligned, "unaligned", false, "Considers sequences as unaligned and format fasta (phylip, nexus,... options are ignored)")
}
