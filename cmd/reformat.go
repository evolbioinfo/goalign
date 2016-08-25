package cmd

import (
	"github.com/spf13/cobra"
	"os"
)

var reformatOutput string
var reformatOutputString string

// reformatCmd represents the reformat command
var reformatCmd = &cobra.Command{
	Use:   "reformat",
	Short: "Reformats input alignment into phylip of fasta format",
	Long: `Reformats input alignment into phylip of fasta format.

Allows to convert formats between phylip and fasta, but also reformat 
fasta and phylip, such as 60 characters per line, etc.
`,
	PersistentPostRun: func(cmd *cobra.Command, args []string) {
		var f *os.File
		var err error

		if reformatOutput == "stdout" || reformatOutput == "-" {
			f = os.Stdout
		} else {
			f, err = os.Create(reformatOutput)
			if err != nil {
				panic(err)
			}
		}
		if rootphylip {
			f.WriteString(reformatOutputString)
		} else {
			f.WriteString(reformatOutputString)
		}
		f.Close()
	},
}

func init() {
	RootCmd.AddCommand(reformatCmd)
	reformatCmd.PersistentFlags().StringVarP(&reformatOutput, "output", "o", "stdout", "Reformated alignment output file")
}
