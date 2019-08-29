package cmd

import (
	"fmt"
	"os"
	"strings"

	"github.com/evolbioinfo/goalign/align"
	"github.com/evolbioinfo/goalign/io"
	"github.com/spf13/cobra"
)

var phaseOutput string
var phaseAAOutput string
var phaseGeneticCode string
var phaseLogOutput string
var orfsequence string
var lencutoff float64
var matchcutoff float64
var phasereverse bool
var phasecutend bool

// translateCmd represents the addid command
var phaseCmd = &cobra.Command{
	Use:   "phase",
	Short: "Find best Starts and set them as new start positions",
	Long: `Find best Starts and set them as new start positions.

This command "phases" input sequences on the basis on either a set of input sequences, or the longest detected orf.
To do so, it will:

1. Search for the longest ORF in the dataset if no reference orf(s) is(are) given;
2. Translate the given ORF(s) in aminoacids;
3. For each sequence of the dataset: translate it in the 3 phases (or 6 if --reverse is given),
   align it with all the translated orfs, and take the phase and the reference orf giving the best alignment;
   If no phase gives a good alignment in any reference orf (cutoffs given by --len-cutoff and --match-cutoff),
   then the sequence flagged as removed;
4. For each sequence, take the Start corresponding to the Start of the ORF, and remove
   nucleotides before (and nucleotides after if --cut-end is given);
5. Return the trimmed nucleotidic sequences (phased), the corresponding amino-acid sequences (phasedaa)
   and the start position on the original nt sequence;
6. The log file contains information on:
    1. Sequence name
    2. Its best matching reference orf
    3. Start position on original nt sequence
    4. Extracted sequence length
    5. Position of the first stop in phase


if --unaligned is set, format options are ignored (phylip, nexus, etc.), and
only Fasta is accepted. Otherwise, alignment is first "unaligned".

If input sequences are not nucleotidic, then returns an error.

Output file is an unaligned set of sequences in fasta.
`,
	RunE: func(cmd *cobra.Command, args []string) (err error) {
		var f, aaf, logf *os.File
		var phased chan align.PhasedSequence
		var inseqs align.SeqBag
		var reforf align.SeqBag
		var orf align.Sequence
		var geneticcode int

		if f, err = openWriteFile(phaseOutput); err != nil {
			io.LogError(err)
			return
		}
		defer closeWriteFile(f, phaseOutput)

		if aaf, err = openWriteFile(phaseAAOutput); err != nil {
			io.LogError(err)
			return
		}
		defer closeWriteFile(aaf, phaseAAOutput)

		if unaligned {
			if inseqs, err = readsequences(infile); err != nil {
				io.LogError(err)
				return
			}
		} else {
			var aligns *align.AlignChannel

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
			if reforf.NbSequences() < 1 {
				err = fmt.Errorf("Reference ORF file should contain at least one sequence")
				io.LogError(err)
				return
			}
		} else {
			// We detect the orf
			if orf, err = inseqs.LongestORF(phasereverse); err != nil {
				io.LogError(err)
				return
			}
			orf.SetName(orf.Name() + "_LongestORF")
			reforf = align.NewSeqBag(align.UNKNOWN)
			reforf.AddSequenceChar(orf.Name(), orf.SequenceChar(), orf.Comment())
			reforf.AutoAlphabet()
		}

		switch phaseGeneticCode {
		case "standard":
			geneticcode = align.GENETIC_CODE_STANDARD
		case "mitov":
			geneticcode = align.GENETIC_CODE_VETEBRATE_MITO
		case "mitoi":
			geneticcode = align.GENETIC_CODE_INVETEBRATE_MITO
		default:
			err = fmt.Errorf("Unknown genetic code : %s", phaseGeneticCode)
			return
		}

		phaser := align.NewPhaser()
		phaser.SetLenCutoff(lencutoff)
		phaser.SetMatchCutoff(matchcutoff)
		phaser.SetReverse(phasereverse)
		phaser.SetCutEnd(phasecutend)
		phaser.SetCpus(rootcpus)
		phaser.SetTranslate(true, geneticcode)
		phaser.SetGapOpen(gapopen)
		phaser.SetGapExtend(gapextend)

		if cmd.Flags().Changed("mismatch") || cmd.Flags().Changed("match") {
			phaser.SetAlignScores(match, mismatch)
		}

		if phased, err = phaser.Phase(reforf, inseqs); err != nil {
			io.LogError(err)
			return
		}

		if logf, err = openWriteFile(phaseLogOutput); err != nil {
			io.LogError(err)
			return
		}
		defer closeWriteFile(logf, phaseLogOutput)

		fmt.Fprintf(logf, "Detected/Given ORF: %s\n", reforf.String())
		fmt.Fprintf(logf, "SeqName\tBestRef\tStartPosition\tExtractedSequenceLength\tFirstStop\n")
		phasedseqs := align.NewSeqBag(align.UNKNOWN)
		phasedseqsaa := align.NewSeqBag(align.UNKNOWN)
		for p := range phased {
			if p.Err != nil {
				err = p.Err
				io.LogError(p.Err)
				return
			}
			if p.Removed {
				fmt.Fprintf(logf, "%s\tN/A\tRemoved\tN/A\n", p.NtSeq.Name())
			} else {
				phasedseqs.AddSequence(p.NtSeq.Name(), p.NtSeq.Sequence(), p.NtSeq.Comment())
				phasedseqsaa.AddSequence(p.AaSeq.Name(), p.AaSeq.Sequence(), p.AaSeq.Comment())
				fmt.Fprintf(logf, "%s\t%s\t%d\t%d\t%d\n", p.NtSeq.Name(), p.Ali.Sequences()[0].Name(), p.Position, p.AaSeq.Length(), strings.Index(p.AaSeq.Sequence(), "*"))
			}
		}

		writeSequences(phasedseqs, f)
		writeSequences(phasedseqsaa, aaf)

		return
	},
}

func init() {
	RootCmd.AddCommand(phaseCmd)
	phaseCmd.PersistentFlags().StringVarP(&phaseOutput, "output", "o", "stdout", "Output \"phased\" FASTA file")
	phaseCmd.PersistentFlags().StringVar(&phaseGeneticCode, "genetic-code", "standard", "Genetic Code: standard, mitoi (invertebrate mitochondrial) or mitov (vertebrate mitochondrial)")
	phaseCmd.PersistentFlags().StringVar(&phaseAAOutput, "aa-output", "none", "Output Met \"phased\" aa FASTA file")
	phaseCmd.PersistentFlags().StringVarP(&phaseLogOutput, "log", "l", "none", "Output log: positions of the considered Start for each sequence")
	phaseCmd.PersistentFlags().Float64Var(&lencutoff, "len-cutoff", -1.0, "Length cutoff, over orf length, to consider sequence hits (-1==No cutoff)")
	phaseCmd.PersistentFlags().Float64Var(&matchcutoff, "match-cutoff", .5, "Nb Matches cutoff, over alignment length, to consider sequence hits (-1==No cutoff)")
	phaseCmd.PersistentFlags().Float64Var(&match, "match", 1.0, "Score for a match for pairwise alignment (if omitted, then take substitution matrix)")
	phaseCmd.PersistentFlags().Float64Var(&mismatch, "mismatch", -1.0, "Score for a mismatch for pairwise alignment (if omitted, then take substitution matrix)")
	phaseCmd.PersistentFlags().Float64Var(&gapopen, "gap-open", -10.0, "Score for opening a gap ")
	phaseCmd.PersistentFlags().Float64Var(&gapextend, "gap-extend", -0.5, "Score for extending a gap ")
	phaseCmd.PersistentFlags().BoolVar(&unaligned, "unaligned", false, "Considers sequences as unaligned and only format fasta is accepted (phylip, nexus,... options are ignored)")
	phaseCmd.PersistentFlags().BoolVar(&phasereverse, "reverse", false, "Search ALSO in the reverse strand (in addition to the forward strand)")
	phaseCmd.PersistentFlags().BoolVar(&phasecutend, "cut-end", false, "Iftrue, then also remove the end of sequences that do not align with orf")
	phaseCmd.PersistentFlags().StringVar(&orfsequence, "ref-orf", "none", "Reference ORF to phase against (if none is given, then will try to get the longest orf in the input data)")
}
