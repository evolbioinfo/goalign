package cmd

import (
	"github.com/evolbioinfo/goalign/align"
	"github.com/evolbioinfo/goalign/io"
	"github.com/evolbioinfo/goalign/io/utils"
	"github.com/spf13/cobra"
)

var sampleseqOutput string
var sampleseqSize int
var sampleseqNbSamples int

// sampleCmd represents the sample command
var sampleseqCmd = &cobra.Command{
	Use:   "seqs",
	Short: "Samples a subset of sequences from the input alignment",
	Long: `Samples a subset of sequences from the input alignment.

May take a Fasta or Phylip alignment in input.

If the input alignment contains several alignments, will process all of them

As output, writes an alignment containing a sample of the sequences.

All alignments are written on the output in the following order:
For each input alignment
  For each sample
    write(sampled alignment)

It is advised to manipulate phylip alignments, in order to be able to
divide the output file with 'goalign divide' for example.

`,
	RunE: func(cmd *cobra.Command, args []string) (err error) {
		var f utils.StringWriterCloser

		if f, err = utils.OpenWriteFile(sampleseqOutput); err != nil {
			io.LogError(err)
			return
		}
		defer utils.CloseWriteFile(f, sampleseqOutput)

		if unaligned {
			var seqs align.SeqBag
			var sample align.SeqBag

			if seqs, err = readsequences(infile); err != nil {
				io.LogError(err)
				return
			}
			for i := 0; i < sampleseqNbSamples; i++ {
				if sample, err = seqs.SampleSeqBag(sampleseqSize); err != nil {
					io.LogError(err)
					return
				}
				writeSequences(sample, f)
			}
		} else {
			var aligns *align.AlignChannel
			var sample align.Alignment

			if aligns, err = readalign(infile); err != nil {
				io.LogError(err)
				return
			}

			for al := range aligns.Achan {
				for i := 0; i < sampleseqNbSamples; i++ {
					if sample, err = al.Sample(sampleseqSize); err != nil {
						io.LogError(err)
						return
					}
					writeAlign(sample, f)
				}
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
	sampleCmd.AddCommand(sampleseqCmd)
	sampleseqCmd.PersistentFlags().BoolVar(&unaligned, "unaligned", false, "Considers sequences as unaligned and format fasta (phylip, nexus,... options are ignored)")
	sampleseqCmd.PersistentFlags().IntVarP(&sampleseqSize, "nb-seq", "n", 1, "Number of sequences to sample from the alignment")
	sampleseqCmd.PersistentFlags().IntVarP(&sampleseqNbSamples, "nb-samples", "s", 1, "Number of samples to generate")
	sampleseqCmd.PersistentFlags().StringVarP(&sampleseqOutput, "output", "o", "stdout", "Sampled alignment output file")
}
