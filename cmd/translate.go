package cmd

import (
	"github.com/fredericlemoine/goalign/io"
	"github.com/spf13/cobra"
)

var translatePhase int
var translateOutput string

// translateCmd represents the addid command
var translateCmd = &cobra.Command{
	Use:   "translate",
	Short: "Translates an input alignment in amino acids",
	Long: `Translates an input alignment in amino acids.

If the input alignment is not nucleotides, then returns an error.

It is possible to drop a given number of characters from the start 
of the alignment, by specifying the '--phase' option.
`,
	Run: func(cmd *cobra.Command, args []string) {
		f := openWriteFile(translateOutput)
		for al := range rootaligns.Achan {
			transAl, err := al.Translate(translatePhase)
			if err != nil {
				io.ExitWithMessage(err)
			}
			writeAlign(transAl, f)
		}
		f.Close()
	},
}

func init() {
	RootCmd.AddCommand(translateCmd)
	translateCmd.PersistentFlags().StringVarP(&translateOutput, "output", "o", "stdout", "Output translated alignment file")
	translateCmd.PersistentFlags().IntVar(&translatePhase, "phase", 0, "Number of characters to drop from the start of the alignment")
}
