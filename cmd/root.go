package cmd

import (
	"bufio"
	"compress/gzip"
	"errors"
	"fmt"
	goio "io"
	"os"
	"runtime"
	"strings"

	"github.com/fredericlemoine/goalign/align"
	"github.com/fredericlemoine/goalign/io"
	"github.com/fredericlemoine/goalign/io/clustal"
	"github.com/fredericlemoine/goalign/io/fasta"
	"github.com/fredericlemoine/goalign/io/nexus"
	"github.com/fredericlemoine/goalign/io/paml"
	"github.com/fredericlemoine/goalign/io/phylip"
	"github.com/fredericlemoine/goalign/io/utils"
	"github.com/fredericlemoine/goalign/version"
	"github.com/spf13/cobra"
)

var cfgFile string
var infile string
var rootaligns align.AlignChannel
var rootphylip bool
var rootnexus bool
var rootclustal bool
var rootcpus int
var rootinputstrict bool = false
var rootoutputstrict bool = false
var rootAutoDetectInputFormat bool

var helptemplate string = `{{with or .Long .Short }}{{. | trim}}

{{end}}Version: ` + version.Version + `

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

Input alignment file formats:
1. Fasta (default)
2. Phylip (-p option)
3. Nexus (-x option)
4. Clustal (-u option)
4. Auto detect (--auto-detect option). In that case, it will test input formats in the following order:
    1. Fasta
    2. Nexus
    3. Phylip
    4. Clustal
    If none of these formats is recognized, then will exit with an error 

Please note that in --auto-detect mode, phylip format is considered as not strict!
`,
	PersistentPreRun: func(cmd *cobra.Command, args []string) {
		runtime.GOMAXPROCS(rootcpus)
		rootaligns = readalign(infile)
	},
	PersistentPostRun: func(cmd *cobra.Command, args []string) {
		if rootaligns.Err != nil {
			io.ExitWithMessage(rootaligns.Err)
		}
	},
}

func readalign(file string) align.AlignChannel {
	var fi goio.Closer
	var r *bufio.Reader
	var err error
	var format int
	var alchan align.AlignChannel

	if fi, r, err = utils.GetReader(file); err != nil {
		io.ExitWithMessage(err)
	}
	if rootAutoDetectInputFormat {
		if alchan, format, err = utils.ParseMultiAlignmentsAuto(fi, r, rootinputstrict); err != nil {
			io.ExitWithMessage(err)
		}
		if format == align.FORMAT_PHYLIP {
			rootphylip = true
		} else if format == align.FORMAT_NEXUS {
			rootnexus = true
		} else if format == align.FORMAT_CLUSTAL {
			rootclustal = true
		}
	} else {
		alchan.Achan = make(chan align.Alignment, 15)
		if rootphylip {
			go func() {
				if err := phylip.NewParser(r, rootinputstrict).ParseMultiple(alchan.Achan); err != nil {
					io.ExitWithMessage(err)
				}
				fi.Close()
			}()
		} else if rootnexus {
			var al align.Alignment
			if al, err = nexus.NewParser(r).Parse(); err != nil {
				io.ExitWithMessage(err)
			}
			alchan.Achan <- al
			fi.Close()
			close(alchan.Achan)
		} else if rootclustal {
			var al align.Alignment
			if al, err = clustal.NewParser(r).Parse(); err != nil {
				io.ExitWithMessage(err)
			}
			alchan.Achan <- al
			fi.Close()
			close(alchan.Achan)
		} else {
			var al align.Alignment
			al, err = fasta.NewParser(r).Parse()
			if err != nil {
				io.ExitWithMessage(err)
			}
			alchan.Achan <- al
			fi.Close()
			close(alchan.Achan)
		}
	}
	return alchan
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
	RootCmd.PersistentFlags().BoolVarP(&rootphylip, "phylip", "p", false, "Alignment is in phylip? default fasta")
	RootCmd.PersistentFlags().BoolVarP(&rootnexus, "nexus", "x", false, "Alignment is in nexus? default fasta")
	RootCmd.PersistentFlags().BoolVarP(&rootclustal, "clustal", "u", false, "Alignment is in clustal? default fasta")

	RootCmd.PersistentFlags().IntVarP(&rootcpus, "threads", "t", 1, "Number of threads")
	RootCmd.PersistentFlags().BoolVar(&rootinputstrict, "input-strict", false, "Strict phylip input format (only used with -p)")
	RootCmd.PersistentFlags().BoolVar(&rootoutputstrict, "output-strict", false, "Strict phylip output format (only used with -p)")

	RootCmd.PersistentFlags().BoolVar(&rootAutoDetectInputFormat, "auto-detect", false, "Auto detects input format (overrides -p, -x and -u)")

	RootCmd.SetHelpTemplate(helptemplate)
}

// initConfig reads in config file and ENV variables if set.
func initConfig() {

}

func writeAlign(al align.Alignment, f *os.File) {
	if rootphylip {
		f.WriteString(phylip.WriteAlignment(al, rootoutputstrict))
	} else if rootnexus {
		f.WriteString(nexus.WriteAlignment(al))
	} else if rootclustal {
		f.WriteString(clustal.WriteAlignment(al))
	} else {
		f.WriteString(fasta.WriteAlignment(al))
	}
}

func writeAlignFasta(al align.Alignment, f *os.File) {
	f.WriteString(fasta.WriteAlignment(al))
}

func writeAlignPhylip(al align.Alignment, f *os.File) {
	f.WriteString(phylip.WriteAlignment(al, rootoutputstrict))
}

func writeAlignNexus(al align.Alignment, f *os.File) {
	f.WriteString(nexus.WriteAlignment(al))
}

func writeAlignClustal(al align.Alignment, f *os.File) {
	f.WriteString(clustal.WriteAlignment(al))
}

func writeUnAlignFasta(al align.Alignment, f *os.File) {
	f.WriteString(fasta.WriteSequences(al))
}

func writeAlignPaml(al align.Alignment, f *os.File) {
	f.WriteString(paml.WriteAlignment(al))
}

func openWriteFile(file string) *os.File {
	if file == "stdout" || file == "-" {
		return os.Stdout
	} else {
		f, err := os.Create(file)
		if err != nil {
			io.ExitWithMessage(err)
		}
		return f
	}
}

// Readln returns a single line (without the ending \n)
// from the input buffered reader.
// An error is returned iff there is an error with the
// buffered reader.
func Readln(r *bufio.Reader) (string, error) {
	var (
		isPrefix bool  = true
		err      error = nil
		line, ln []byte
	)
	for isPrefix && err == nil {
		line, isPrefix, err = r.ReadLine()
		ln = append(ln, line...)
	}
	return string(ln), err
}

func readMapFile(file string, revert bool) (map[string]string, error) {
	outmap := make(map[string]string, 0)
	var mapfile *os.File
	var err error
	var reader *bufio.Reader

	if mapfile, err = os.Open(file); err != nil {
		return outmap, err
	}

	if strings.HasSuffix(file, ".gz") {
		if gr, err2 := gzip.NewReader(mapfile); err2 != nil {
			return outmap, err2
		} else {
			reader = bufio.NewReader(gr)
		}
	} else {
		reader = bufio.NewReader(mapfile)
	}
	line, e := Readln(reader)
	nl := 1
	for e == nil {
		cols := strings.Split(line, "\t")
		if len(cols) != 2 {
			return outmap, errors.New("Map file does not have 2 fields at line: " + fmt.Sprintf("%d", nl))
		}
		if revert {
			outmap[cols[1]] = cols[0]
		} else {
			outmap[cols[0]] = cols[1]
		}
		line, e = Readln(reader)
		nl++
	}

	if err = mapfile.Close(); err != nil {
		return outmap, err
	}

	return outmap, nil
}
