package cmd

import (
	"fmt"
	"os"

	"github.com/evolbioinfo/goalign/align"
	"github.com/evolbioinfo/goalign/io"
	"github.com/spf13/cobra"
)

var divideOutput string
var divideoutputFasta bool
var divideNbSeqs int
var divideUnaligned bool

// divideCmd represents the divide command
var divideCmd = &cobra.Command{
	Use:   "divide",
	Short: "Divide an input alignment in several output files",
	Long: `Divide an input alignment in several output files

The default behavior is to take an input alignment file containing 
potentially several alignments (e.g. with Phylip format ), and output
one alignment per output file.

If the alignment is in fasta format : will create 1 file
Otherwise, will create one output file per input alignment.

if the option --nb-sequences <n> is given, then outputs n sequences
per output file. 

-o : is the prefix of output files
if -o div, it will create files div_0.ph...div_n.ph

Output files will be in Phylip Format or in fasta format depending on -f

Example:

gotree divide -i align.ph -p -o out

`,
	RunE: func(cmd *cobra.Command, args []string) (err error) {
		var aligns *align.AlignChannel
		var tmpAlign align.Alignment
		var tmpSeqs align.SeqBag

		var f *os.File

		i := 0

		ext := alignExtension()
		if divideoutputFasta {
			ext = ".fa"
		}

		if divideUnaligned {
			var seqs align.SeqBag
			if seqs, err = readsequences(infile); err != nil {
				io.LogError(err)
				return
			}

			if divideNbSeqs == 0 {
				if f, err = openWriteFile(fmt.Sprintf("%s_%03d%s", divideOutput, i, ext)); err != nil {
					io.LogError(err)
					return
				}
				writeSequences(seqs, f)
				f.Close()
				i++
			} else {
				tmpSeqs = align.NewSeqBag(seqs.Alphabet())
				nb := 0
				seqs.IterateAll(func(name string, sequence []uint8, comment string) bool {
					tmpSeqs.AddSequenceChar(name, sequence, comment)
					nb++
					if divideNbSeqs > 0 && nb%divideNbSeqs == 0 {
						if f, err = openWriteFile(fmt.Sprintf("%s_%03d%s", divideOutput, i, ext)); err != nil {
							io.LogError(err)
							return true
						}
						writeSequences(tmpSeqs, f)
						f.Close()
						tmpSeqs = align.NewSeqBag(seqs.Alphabet())
						i++
					}
					return false
				})
				if nb%divideNbSeqs > 0 {
					if f, err = openWriteFile(fmt.Sprintf("%s_%03d%s", divideOutput, i, ext)); err != nil {
						io.LogError(err)
					}
					writeSequences(tmpSeqs, f)
					f.Close()
					i++
				}
			}
		} else {
			if aligns, err = readalign(infile); err != nil {
				io.LogError(err)
				return
			}

			for al := range aligns.Achan {
				if divideNbSeqs == 0 {
					if f, err = openWriteFile(fmt.Sprintf("%s_%03d%s", divideOutput, i, ext)); err != nil {
						io.LogError(err)
						return
					}
					if divideoutputFasta {
						writeAlignFasta(al, f)
					} else {
						writeAlign(al, f)
					}
					f.Close()
					i++
				} else {
					tmpAlign = align.NewAlign(al.Alphabet())
					nb := 0
					al.IterateAll(func(name string, sequence []uint8, comment string) bool {
						tmpAlign.AddSequenceChar(name, sequence, comment)
						nb++
						if nb%divideNbSeqs == 0 {
							if f, err = openWriteFile(fmt.Sprintf("%s_%03d%s", divideOutput, i, ext)); err != nil {
								io.LogError(err)
								return true
							}
							if divideoutputFasta {
								writeAlignFasta(tmpAlign, f)
							} else {
								writeAlign(tmpAlign, f)
							}
							f.Close()
							i++
							tmpAlign = align.NewAlign(al.Alphabet())
						}
						return false
					})
					if nb%divideNbSeqs > 0 {
						if f, err = openWriteFile(fmt.Sprintf("%s_%03d%s", divideOutput, i, ext)); err != nil {
							io.LogError(err)
						}
						if divideoutputFasta {
							writeAlignFasta(tmpAlign, f)
						} else {
							writeAlign(tmpAlign, f)
						}
						f.Close()
						i++
					}
				}
				f.Close()
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
	RootCmd.AddCommand(divideCmd)
	divideCmd.PersistentFlags().StringVarP(&divideOutput, "output", "o", "prefix", "Divided alignment output files prefix")
	divideCmd.PersistentFlags().IntVar(&divideNbSeqs, "nb-sequences", 0, "Number of sequences per output file (<=0 : all sequences, >0: each alignment will be additionnaly split in several alignments)")
	divideCmd.PersistentFlags().BoolVar(&divideUnaligned, "unaligned", false, "Considers sequences as unaligned and format fasta (phylip, nexus,... options are ignored)")
	divideCmd.PersistentFlags().BoolVarP(&divideoutputFasta, "out-fasta", "f", false, "Forces output files to be in fasta format")

}
