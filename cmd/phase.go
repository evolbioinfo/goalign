package cmd

import (
	"fmt"
	"os"

	"github.com/fredericlemoine/goalign/align"
	"github.com/fredericlemoine/goalign/io"
	"github.com/spf13/cobra"
)

var phaseOutput string
var phaseLogOutput string

// translateCmd represents the addid command
var phaseCmd = &cobra.Command{
	Use:   "phase",
	Short: "Find best ATG and set it as new start position",
	Long: `Find best ATG and set it as new start position.

if --unaligned is set, format options are ingored (phylip, nexus, etc.), and
only Fasta is accepted.

If input sequences are not nucleotidic, then returns an error.

Output file is an unaligned set of sequences in fasta.

`,
	RunE: func(cmd *cobra.Command, args []string) (err error) {
		var f, log *os.File
		var pos []int
		var phased align.SeqBag
		var inseqs align.SeqBag

		if f, err = openWriteFile(phaseOutput); err != nil {
			io.LogError(err)
			return
		}
		defer closeWriteFile(f, phaseOutput)

		if log, err = openWriteFile(phaseLogOutput); err != nil {
			io.LogError(err)
			return
		}
		defer closeWriteFile(log, phaseLogOutput)

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

		if phased, pos, err = inseqs.Phase(); err != nil {
			io.LogError(err)
			return
		}

		if phaseLogOutput != "none" {
			for i, v := range pos {
				n, _ := phased.GetSequenceNameById(i)
				fmt.Fprintf(log, "%s\t%d\n", n, v)
			}
		}

		writeSequences(phased, f)

		return
	},
}

func init() {
	RootCmd.AddCommand(phaseCmd)
	phaseCmd.PersistentFlags().StringVarP(&phaseOutput, "output", "o", "stdout", "Output translated alignment file")
	phaseCmd.PersistentFlags().StringVarP(&phaseLogOutput, "log", "l", "none", "Output log: positions of the considered ATG for each sequence")
	phaseCmd.PersistentFlags().BoolVar(&unaligned, "unaligned", false, "Considers sequences as unaligned and only format fasta is accepted (phylip, nexus,... options are ignored)")
}
