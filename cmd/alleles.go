package cmd

import (
	"fmt"

	"github.com/spf13/cobra"
)

// allelesCmd represents the alleles command
var allelesCmd = &cobra.Command{
	Use:   "alleles",
	Short: "Prints the average number of alleles per sites of the alignment",
	Long:  `Prints the average number of alleles per sites of the alignment.`,
	Run: func(cmd *cobra.Command, args []string) {
		aligns := readalign(infile)
		for al := range aligns.Achan {
			fmt.Println(al.AvgAllelesPerSite())
		}
	},
}

func init() {
	statsCmd.AddCommand(allelesCmd)
}
