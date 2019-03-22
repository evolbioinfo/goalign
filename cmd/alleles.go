package cmd

import (
	"fmt"

	"github.com/evolbioinfo/goalign/align"
	"github.com/evolbioinfo/goalign/io"
	"github.com/spf13/cobra"
)

// allelesCmd represents the alleles command
var allelesCmd = &cobra.Command{
	Use:   "alleles",
	Short: "Prints the average number of alleles per sites of the alignment",
	Long:  `Prints the average number of alleles per sites of the alignment.`,
	RunE: func(cmd *cobra.Command, args []string) (err error) {
		var aligns *align.AlignChannel

		if aligns, err = readalign(infile); err != nil {
			io.LogError(err)
			return
		}
		for al := range aligns.Achan {
			fmt.Println(al.AvgAllelesPerSite())
		}

		if aligns.Err != nil {
			err = aligns.Err
			io.LogError(err)
		}
		return
	},
}

func init() {
	statsCmd.AddCommand(allelesCmd)
}
