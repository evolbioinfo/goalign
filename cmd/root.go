package cmd

import (
	"bufio"
	"compress/gzip"
	"fmt"
	"github.com/fredericlemoine/goalign/align"
	"github.com/fredericlemoine/goalign/io/fasta"
	"github.com/fredericlemoine/goalign/io/phylip"
	"github.com/spf13/cobra"
	"os"
	"strings"
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
		var r *bufio.Reader
		var err error
		if infile == "stdin" || infile == "-" {
			fi = os.Stdin
		} else {
			fi, err = os.Open(infile)
			if err != nil {
				panic(err)
			}
		}
		if strings.HasSuffix(infile, ".gz") {
			if gr, err := gzip.NewReader(fi); err != nil {
				panic(err)
			} else {
				r = bufio.NewReader(gr)
			}

		} else {
			r = bufio.NewReader(fi)
		}
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

func writeAlign(al align.Alignment, file string) {
	var f *os.File
	var err error

	if file == "stdout" || file == "-" {
		f = os.Stdout
	} else {
		f, err = os.Create(file)
		if err != nil {
			panic(err)
		}
	}
	if rootphylip {
		f.WriteString(phylip.WriteAlignment(al))
	} else {
		f.WriteString(fasta.WriteAlignment(al))
	}
	f.Close()
}
