package cmd

import (
	"github.com/spf13/cobra"
)

var mutateOutput string
var mutateRate float64

// mutateCmd represents the mutate command
var mutateCmd = &cobra.Command{
	Use:   "mutate",
	Short: "Adds substitutions (~sequencing errors), or gaps, uniformly in an input alignment",
	Long: `Adds substitutions (~sequencing error), or gaps, uniformly in an input alignment.
`,
}

func init() {
	RootCmd.AddCommand(mutateCmd)
	mutateCmd.PersistentFlags().Float64VarP(&mutateRate, "rate", "r", 0.1, "Mutation rate per nucleotide/amino acid")
	mutateCmd.PersistentFlags().StringVarP(&mutateOutput, "output", "o", "stdout", "Mutated alignment output file")
}
