package cmd

import (
	"fmt"

	"github.com/evolbioinfo/goalign/align"
	"github.com/evolbioinfo/goalign/io"
	"github.com/evolbioinfo/goalign/io/utils"
	"github.com/spf13/cobra"
)

var translatePhase int
var translateOutput string
var translateGeneticCode string
var translaterefseq string

// translateCmd represents the addid command
var translateCmd = &cobra.Command{
	Use:   "translate",
	Short: "Translates an input alignment in amino acids",
	Long: `Translates an input alignment in amino acids.

If the input alignment is not nucleotides, then returns an error.

It is possible to drop a given number of characters from the start 
of the alignment, by specifying the '--phase' option.

If given phase is -1, then it will translate in the 3 phases, 
from positions 0, 1 and 2. Sequence names will be added the suffix
_<phase>. At the end, 3x times more sequences will be present in the
file.

It is possible to specify alternative genetic code with --genetic-code 
(mitoi, mitov, or standard).

IUPAC codes are taken into account for the translation. If a codon containing 
IUPAC code is ambiguous for translation, then a X is added in place of the aminoacid.

If --ref-seq is given, be careful about the behavior! As with goalign extract, it will will translate
the alignment with the following process: The alignment will be translated codon by
codon using the given reference sequence as guide, by iterating over the reference non gap nucleotides 3 by 3. 
At each iteration, the current reference codon may have gaps between nucleotides, and the translation of the
current codon will be done as following:
	* ex 1:
		Ref: AC--GTACGT
		Seq: ACTTGTACGT
		In that case, the first ref codon is [0,1,4], corresponding to sequence ACTTG in seq
		ACTTG % 3 != 0 ==> Frameshift? => Replaced by T in ref and X in the compared sequence.
	* ex 2:
		Ref: AC---GTACGT
		Seq: ACTTTGTACGT
		ref codon: [0,1,5]
		seq      : ACTTTG (%3==0): Insertion - OK => Replaced by "T-" in ref and "TL" in seq
	* ex 3:
		Ref: ACGTACGT
		Seq: A--TACGT
		ref codon: [0,1,2]
		seq      : A--: Deletion: not ok : Frameshift? => Replaced by "T" in ref and "X" in comp
	* ex 4:
		Ref: AC----GTACGT
		Seq: ACTT-TGTACGT
		ref codon: [0,1,6]
		seq      : ACTTTG (%3==0): Insertion - OK => Replaced by "T-" in ref and "TT" in seq
	* ex 5:
		Ref: AC----GTACGT
		Seq: ACT--TGTACGT
		ref codon: [0,1,6]
		seq      : ACTTTG : Insertion not OK : Frameshift? => Replaced by "T-" in ref and "XX" in seq
This allows to easily translate a multiple sequence alignment containing partial sequences, but the 
interpretation should be careful: the translation of some sequences may not be representative of the 
translation of the unaligned sequences.
`,
	RunE: func(cmd *cobra.Command, args []string) (err error) {
		var f utils.StringWriterCloser
		var geneticcode int

		if f, err = utils.OpenWriteFile(translateOutput); err != nil {
			io.LogError(err)
			return
		}
		defer utils.CloseWriteFile(f, translateOutput)

		switch translateGeneticCode {
		case "standard":
			geneticcode = align.GENETIC_CODE_STANDARD
		case "mitov":
			geneticcode = align.GENETIC_CODE_VETEBRATE_MITO
		case "mitoi":
			geneticcode = align.GENETIC_CODE_INVETEBRATE_MITO

		default:
			err = fmt.Errorf("unknown genetic code : %s", translateGeneticCode)
			return
		}

		if unaligned {
			var seqs align.SeqBag

			if seqs, err = readsequences(infile); err != nil {
				io.LogError(err)
				return
			}
			if err = seqs.Translate(translatePhase, geneticcode); err != nil {
				io.LogError(err)
				return
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
				if translaterefseq != "" {
					if err = al.TranslateByReference(translatePhase, geneticcode, translaterefseq); err != nil {
						io.LogError(err)
						return
					}
				} else {
					if err = al.Translate(translatePhase, geneticcode); err != nil {
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
	RootCmd.AddCommand(translateCmd)
	translateCmd.PersistentFlags().StringVar(&translateGeneticCode, "genetic-code", "standard", "Genetic Code: standard, mitoi (invertebrate mitochondrial) or mitov (vertebrate mitochondrial)")
	translateCmd.PersistentFlags().StringVarP(&translateOutput, "output", "o", "stdout", "Output translated alignment file")
	translateCmd.PersistentFlags().StringVar(&translaterefseq, "ref-seq", "", "Reference sequence on which coordinates are given (ignored if --unaligned)")
	translateCmd.PersistentFlags().IntVar(&translatePhase, "phase", 0, "Number of characters to drop from the start of the alignment (if -1: Translate in the 3 phases, from positions 0, 1, and 2)")
	translateCmd.PersistentFlags().BoolVar(&unaligned, "unaligned", false, "Considers sequences as unaligned and format fasta (phylip, nexus,... options are ignored)")
}
