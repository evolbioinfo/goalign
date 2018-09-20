package cmd

import (
	"github.com/fredericlemoine/goalign/io"
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
	Run: func(cmd *cobra.Command, args []string) {
		aligns := readalign(infile)
		f := openWriteFile(dedupOutput)
		for al := range aligns.Achan {
			if out, err := al.Deduplicate(); err != nil {
				if err != nil {
					io.ExitWithMessage(err)
				}
			} else {
				writeAlign(out, f)
			}
		}
		f.Close()
	},
}

func init() {
	RootCmd.AddCommand(dedupCmd)
	dedupCmd.PersistentFlags().StringVarP(&dedupOutput, "output", "o", "stdout", "Deduplicated output alignment file")
}
