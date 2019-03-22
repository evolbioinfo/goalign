package cmd

import (
	"os"

	"github.com/evolbioinfo/goalign/align"
	"github.com/evolbioinfo/goalign/io"
	"github.com/spf13/cobra"
)

var dedupOutput string

// dedupCmd represents the dedup command
var dedupCmd = &cobra.Command{
	Use:   "dedup",
	Short: "Deduplicate sequences that have the same sequence",
	Long: `Deduplicate sequences that have the same sequence

The name of the first sequence is kept

Example: 

ali.phy
1 AAAAAA
2 CCCCCC
3 GGGGGG
4 GGGGGG

goalign dedup -i ali.phy will produce:

1 AAAAAA
2 CCCCCC
3 GGGGGG
`,
	RunE: func(cmd *cobra.Command, args []string) (err error) {
		var aligns *align.AlignChannel
		var f *os.File

		if aligns, err = readalign(infile); err != nil {
			io.LogError(err)
			return
		}
		if f, err = openWriteFile(dedupOutput); err != nil {
			io.LogError(err)
			return
		}
		defer closeWriteFile(f, dedupOutput)

		for al := range aligns.Achan {
			if err = al.Deduplicate(); err != nil {
				io.LogError(err)
				return
			} else {
				writeAlign(al, f)
			}
		}

		if aligns.Err != nil {
			err = aligns.Err
			io.LogError(err)
		}
		return
	},
}

func init() {
	RootCmd.AddCommand(dedupCmd)
	dedupCmd.PersistentFlags().StringVarP(&dedupOutput, "output", "o", "stdout", "Deduplicated output alignment file")
}
