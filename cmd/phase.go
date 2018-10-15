package cmd

import (
	"fmt"
	"log"
	"os"

	"github.com/fredericlemoine/goalign/align"
	"github.com/fredericlemoine/goalign/io"
	"github.com/spf13/cobra"
)

var phaseOutput string
var phaseLogOutput string
var orfsequence string

// translateCmd represents the addid command
var phaseCmd = &cobra.Command{
	Use:   "phase",
	Short: "Find best ATGs and set them as new start positions",
	Long: `Find best ATGs and set them as new start positions.

if --unaligned is set, format options are ignored (phylip, nexus, etc.), and
only Fasta is accepted.

If input sequences are not nucleotidic, then returns an error.

Output file is an unaligned set of sequences in fasta.

`,
	RunE: func(cmd *cobra.Command, args []string) (err error) {
		var f, logf *os.File
		var pos []int
		var phased align.SeqBag
		var inseqs align.SeqBag
		var reforf align.SeqBag
		var orf align.Sequence
		var ok bool

		if f, err = openWriteFile(phaseOutput); err != nil {
			io.LogError(err)
			return
		}
		defer closeWriteFile(f, phaseOutput)

		if logf, err = openWriteFile(phaseLogOutput); err != nil {
			io.LogError(err)
			return
		}
		defer closeWriteFile(logf, phaseLogOutput)

		if orfsequence != "none" {
			if reforf, err = readsequences(orfsequence); err != nil {
				io.LogError(err)
				return
			}
			if reforf.NbSequences() != 1 {
				err = fmt.Errorf("Reference ORF file should contain only one sequence")
				io.LogError(err)
				return
			}
			if orf, ok = reforf.Sequence(0); !ok {
				io.LogError(fmt.Errorf("Sequence 0 is not present in the orf file"))
				return
			}
		}

		if unaligned {
			if inseqs, err = readsequences(infile); err != nil {
				io.LogError(err)
				return
			}
		} else {
			var aligns align.AlignChannel

			if aligns, err = readalign(infile); err != nil {
				io.LogError(err)
				return
			}
			inseqs = (<-aligns.Achan).Unalign()
		}

		if phased, pos, err = inseqs.Phase(orf); err != nil {
			io.LogError(err)
			return
		}

		if phaseLogOutput != "none" {
			for i, v := range pos {
				n, _ := phased.GetSequenceNameById(i)
				fmt.Fprintf(logf, "%s\t%d\n", n, v)
			}
		}

		writeSequences(phased, f)

		log.SetOutput(os.Stderr)

		return
	},
}

func init() {
	RootCmd.AddCommand(phaseCmd)
	phaseCmd.PersistentFlags().StringVarP(&phaseOutput, "output", "o", "stdout", "Output ATG \"phased\" FASTA file")
	phaseCmd.PersistentFlags().StringVarP(&phaseLogOutput, "log", "l", "none", "Output log: positions of the considered ATG for each sequence")
	phaseCmd.PersistentFlags().BoolVar(&unaligned, "unaligned", false, "Considers sequences as unaligned and only format fasta is accepted (phylip, nexus,... options are ignored)")
	phaseCmd.PersistentFlags().StringVar(&orfsequence, "ref-orf", "none", "Reference ORF to phase against (if none is given, then will try to get the longest orf in the input data)")
}
