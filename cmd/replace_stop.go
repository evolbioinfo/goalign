package cmd

import (
	"fmt"

	"github.com/evolbioinfo/goalign/align"
	"github.com/evolbioinfo/goalign/io"
	"github.com/evolbioinfo/goalign/io/utils"
	"github.com/spf13/cobra"
)

var replacePhase int
var replaceGeneticCode string

// renameCmd represents the rename command
var replaceStopCmd = &cobra.Command{
	Use:   "stops",
	Short: "Replace STOP codons in input nt sequences by NNN",
	Long: `Replace STOP codons in input nt sequences by NNN, in the given phase (except the last codon).

	If --phase is given (>=0), then starts at the given offset (default=0)
`,
	RunE: func(cmd *cobra.Command, args []string) (err error) {
		var aligns *align.AlignChannel
		var f utils.StringWriterCloser
		var seqs align.SeqBag
		var geneticcode int

		if replacePhase < 0 {
			err = fmt.Errorf("phase must be >=0")
			return
		}

		switch replaceGeneticCode {
		case "standard":
			geneticcode = align.GENETIC_CODE_STANDARD
		case "mitov":
			geneticcode = align.GENETIC_CODE_VETEBRATE_MITO
		case "mitoi":
			geneticcode = align.GENETIC_CODE_INVETEBRATE_MITO

		default:
			err = fmt.Errorf("unknown genetic code : %s", replaceGeneticCode)
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
			if err = seqs.ReplaceStops(replacePhase, geneticcode); err != nil {
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
				if err = al.ReplaceStops(replacePhase, geneticcode); err != nil {
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
	replaceCmd.AddCommand(replaceStopCmd)

	replaceStopCmd.PersistentFlags().IntVar(&replacePhase, "phase", 0, "Phase in which replace STOP codons")
	replaceStopCmd.PersistentFlags().StringVar(&replaceGeneticCode, "genetic-code", "standard", "Genetic Code: standard, mitoi (invertebrate mitochondrial) or mitov (vertebrate mitochondrial)")

}
