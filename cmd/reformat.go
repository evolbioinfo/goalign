package cmd

import (
	"github.com/spf13/cobra"
)

var reformatOutput string
var reformatCleanNames bool

// reformatCmd represents the reformat command
var reformatCmd = &cobra.Command{
	Use:   "reformat",
	Short: "Reformats input alignment into phylip of fasta format",
	Long: `Reformats input alignment into phylip of fasta format.

Allows to convert formats between phylip, fasta and nexus, but also reformats
fasta and phylip, such as 60 characters per line, etc.

`,
}

func init() {
	RootCmd.AddCommand(reformatCmd)
	reformatCmd.PersistentFlags().StringVarP(&reformatOutput, "output", "o", "stdout", "Reformated alignment output file")
	reformatCmd.PersistentFlags().BoolVar(&reformatCleanNames, "clean-names", false, "Replaces special characters (tabs, spaces, newick characters) with '-' from input sequence names before writing output alignment")
}
