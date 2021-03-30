package cmd

import (
	"os"

	"github.com/evolbioinfo/goalign/align"
	"github.com/evolbioinfo/goalign/io"
	"github.com/spf13/cobra"
)

var consensusOutput string
var consensusExcludeGaps bool
var consensusIgnoreGaps bool
var consensusIgnoreNs bool

// concatCmd represents the concat command
var consensusCmd = &cobra.Command{
	Use:   "consensus",
	Short: "Compute the majority consensus of an input alignment",
	Long: `Compute the majority consensus of an input alignment.

For example:

goalign consensus -i align.phylip -p 

It will generate a single sequence whose sites will correspond to the
majority characters at each positions (including gaps).

If several alignment are present in the input file (for phylip)
then will output several consensus sequences.
`,
	RunE: func(cmd *cobra.Command, args []string) (err error) {
		var aligns *align.AlignChannel
		var f *os.File

		consensusIgnoreGaps = consensusIgnoreGaps || consensusExcludeGaps

		if aligns, err = readalign(infile); err != nil {
			io.LogError(err)
			return
		}

		if f, err = openWriteFile(consensusOutput); err != nil {
			io.LogError(err)
			return
		}
		defer closeWriteFile(f, consensusOutput)

		for al := range aligns.Achan {
			cons := al.Consensus(consensusIgnoreGaps, consensusIgnoreNs)
			writeAlign(cons, f)
		}

		if aligns.Err != nil {
			err = aligns.Err
			io.LogError(err)
		}
		return
	},
}

func init() {
	RootCmd.AddCommand(consensusCmd)
	consensusCmd.PersistentFlags().StringVarP(&consensusOutput, "output", "o", "stdout", "Alignment output file")
	consensusCmd.PersistentFlags().BoolVar(&consensusIgnoreGaps, "ignore-gaps", false, "Ignore gaps in the majority computation")
	consensusCmd.PersistentFlags().BoolVar(&consensusExcludeGaps, "exclude-gaps", false, "Ignore gaps in the majority computation (for backward compatibility, will be removed in future releases)")
	consensusCmd.PersistentFlags().BoolVar(&consensusIgnoreNs, "ignore-n", false, "Ignore Ns in the majority computation")
}
