package cmd

import (
	"bufio"
	"fmt"
	"github.com/fredericlemoine/goalign/align"
	"github.com/fredericlemoine/goalign/io/fasta"
	"github.com/fredericlemoine/goalign/io/phylip"
	"github.com/spf13/cobra"
	"os"
)

var cfgFile string
var infile string
var rootalign align.Alignment
var rootphylip bool

var Version string = "Unset"

var helptemplate string = `{{with or .Long .Short }}{{. | trim}}

{{end}}Version: ` + Version + `

{{if or .Runnable .HasSubCommands}}{{.UsageString}}{{end}}
`

// RootCmd represents the base command when called without any subcommands
var RootCmd = &cobra.Command{
	Use:   "goalign",
	Short: "goalign: A set of tools to handle sequence alignments",
	Long: `goalign: A set of tools to handle sequence alignments.

It allows to :
1 - Convert formats
2 - Shuffle alignments
3 - Sample alignments
4 - Print informations about alignments
5 - Generate bootstrap alignments

`,
	PersistentPreRun: func(cmd *cobra.Command, args []string) {
		var fi *os.File
		var err error
		if infile == "stdin" {
			fi = os.Stdin
		} else {
			fi, err = os.Open(infile)
			if err != nil {
				panic(err)
			}
		}
		r := bufio.NewReader(fi)
		var al align.Alignment
		var err2 error
		if rootphylip {
			al, err2 = phylip.NewParser(r).Parse()
		} else {
			al, err2 = fasta.NewParser(r).Parse()
		}
		if err2 != nil {
			panic(err2)
		}
		rootalign = al
	},
}

// Execute adds all child commands to the root command sets flags appropriately.
// This is called by main.main(). It only needs to happen once to the rootCmd.
func Execute() {
	if err := RootCmd.Execute(); err != nil {
		fmt.Println(err)
		os.Exit(-1)
	}
}

func init() {
	cobra.OnInitialize(initConfig)

	RootCmd.PersistentFlags().StringVarP(&infile, "align", "i", "stdin", "Alignment input file")
	RootCmd.PersistentFlags().BoolVarP(&rootphylip, "phylip", "p", false, "Alignment is in phylip? False=Fasta")
	RootCmd.SetHelpTemplate(helptemplate)
}

// initConfig reads in config file and ENV variables if set.
func initConfig() {

}
