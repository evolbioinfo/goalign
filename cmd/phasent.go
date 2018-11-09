package cmd

import (
	"fmt"
	"os"

	"github.com/fredericlemoine/goalign/align"
	"github.com/fredericlemoine/goalign/io"
	"github.com/spf13/cobra"
)

// translateCmd represents the addid command
var phasentCmd = &cobra.Command{
	Use:   "phasent",
	Short: "Find best ATGs and set them as new start positions",
	Long: `Find best ATGs and set them as new start positions.

if --unaligned is set, format options are ignored (phylip, nexus, etc.), and
only Fasta is accepted.

Unlike goalign phase, it does not take into account protein information (3 or 6 phases).

If input sequences are not nucleotidic, then returns an error.

Output files:

--output : unaligned set of phased sequences in fasta


 1. If alignment is bad (>lencutoff * orf length, >matchcutoff matches over the 
    align length and starting at first position of the ORF), then the sequence 
    is discarded;
 3. For each sequence, take the Start corresponding to the Start of the ORF, and
    remove nucleotides before;
 4. Return the trimmed nucleotidic sequences (phased), the positions of starts in
    the nucleotidic sequences, and the removed sequence names.
 If cutend is true, then also remove the end of sequences that do not align with orf
 It does not modify the input object
`,
	RunE: func(cmd *cobra.Command, args []string) (err error) {
		var f, logf *os.File
		var pos []int
		var phased align.SeqBag
		var inseqs align.SeqBag
		var reforf align.SeqBag
		var orf align.Sequence
		var ok bool
		var removed []string

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
		} else {
			// We detect the orf
			if orf, err = inseqs.LongestORF(phasereverse); err != nil {
				io.LogError(err)
				return
			}
			orf.SetName(orf.Name() + "_LongestORF")
		}

		if phased, pos, removed, err = inseqs.PhaseNt(orf, lencutoff, matchcutoff, phasereverse, phasecutend, rootcpus); err != nil {
			io.LogError(err)
			return
		}

		if phaseLogOutput != "none" {
			fmt.Fprintf(logf, "Detected/Given ORF in %s: %s\n", orf.Name(), orf.Sequence())
			for i, v := range pos {
				n, _ := phased.GetSequenceNameById(i)
				fmt.Fprintf(logf, "%s\t%d\n", n, v)
			}
			for _, v := range removed {
				fmt.Fprintf(logf, "Removed: %s\n", v)
			}

		}

		writeSequences(phased, f)

		return
	},
}

func init() {
	RootCmd.AddCommand(phasentCmd)
	phasentCmd.PersistentFlags().StringVarP(&phaseOutput, "output", "o", "stdout", "Output ATG \"phased\" FASTA file")
	phasentCmd.PersistentFlags().StringVarP(&phaseLogOutput, "log", "l", "none", "Output log: positions of the considered ATG for each sequence")
	phasentCmd.PersistentFlags().Float64Var(&lencutoff, "len-cutoff", -1.0, "Length cutoff, over orf length, to consider sequence hits (-1==No cutoff)")
	phasentCmd.PersistentFlags().Float64Var(&matchcutoff, "match-cutoff", .5, "Nb Matches cutoff, over alignment length, to consider sequence hits (-1==No cutoff)")
	phasentCmd.PersistentFlags().BoolVar(&unaligned, "unaligned", false, "Considers sequences as unaligned and only format fasta is accepted (phylip, nexus,... options are ignored)")
	phasentCmd.PersistentFlags().BoolVar(&phasereverse, "reverse", false, "Search ALSO in the reverse strand (in addition to the forward strand)")
	phasentCmd.PersistentFlags().BoolVar(&phasecutend, "cut-end", false, "Iftrue, then also remove the end of sequences that do not align with orf")
	phasentCmd.PersistentFlags().StringVar(&orfsequence, "ref-orf", "none", "Reference ORF to phase against (if none is given, then will try to get the longest orf in the input data)")
}
