package cmd

import (
	"os"

	"github.com/fredericlemoine/goalign/align"
	"github.com/fredericlemoine/goalign/io"
	"github.com/spf13/cobra"
)

var siteRate float64
var siteRogue float64
var siteRogueNameFile string
var stableRogues bool

// sitesCmd represents the sites command
var sitesCmd = &cobra.Command{
	Use:   "sites",
	Short: "Shuffles n alignment sites vertically",
	Long: `Shuffles n alignment sites vertically.

It may take Fasta or Phylip input alignment.

If the input alignment contains several alignments, will process all of them

The only option is to specify the rate of shuffled sites.
A rate of 0.5 will shuffle 50% of the sites of the alignment.

It takes n sites of the input alignment and reassign the 
characters to different sequences.

Example of usage:

goalign shuffle sites -i align.phylip -p -r 0.5
goalign shuffle sites -i align.fasta -r 0.5

`,
	RunE: func(cmd *cobra.Command, args []string) (err error) {
		var aligns align.AlignChannel
		var f *os.File
		var nameFile *os.File

		aligns, err = readalign(infile)
		if f, err = openWriteFile(shuffleOutput); err != nil {
			io.LogError(err)
			return
		}
		defer closeWriteFile(f, shuffleOutput)

		if nameFile, err = openWriteFile(siteRogueNameFile); err != nil {
			io.LogError(err)
			return
		}
		defer closeWriteFile(nameFile, siteRogueNameFile)

		for al := range aligns.Achan {
			names := al.ShuffleSites(siteRate, siteRogue, stableRogues)
			writeAlign(al, f)
			for _, n := range names {
				nameFile.WriteString(n)
				nameFile.WriteString("\n")
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
	shuffleCmd.AddCommand(sitesCmd)
	sitesCmd.PersistentFlags().Float64VarP(&siteRate, "rate", "r", 0.5, "Rate of shuffled sites (>=0 and <=1)")
	sitesCmd.PersistentFlags().Float64Var(&siteRogue, "rogue", 0.0, "If set, then will take the given proportion of taxa, and will apply shuffle again on --rate of the remaining intact sites")
	sitesCmd.PersistentFlags().StringVar(&siteRogueNameFile, "rogue-file", "stdout", "Rogue sequence names output file")
	sitesCmd.PersistentFlags().BoolVar(&stableRogues, "stable-rogues", false, "If true, then with a given seed, rogues will always be the same with all alignments having sequences in the same order. It may not be the case if false, especially when alignemnts have different lengths.")
}
