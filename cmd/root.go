package cmd

import (
	"bufio"
	"compress/gzip"
	"errors"
	"fmt"
	goio "io"
	"math/rand"
	"os"
	"runtime"
	"strings"
	"time"

	"github.com/fredericlemoine/cobrashell"
	"github.com/fredericlemoine/goalign/align"
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
var rootphylip bool
var rootnexus bool
var rootclustal bool
var rootcpus int
var rootinputstrict bool = false
var rootoutputstrict bool = false
var rootoutputoneline = false
var rootoutputnoblock = false
var rootAutoDetectInputFormat bool
var seed int64 = -1
var unaligned bool

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
		if seed == -1 {
			seed = time.Now().UTC().UnixNano()
		}
		rand.Seed(seed)
	},
	Run: func(cmd *cobra.Command, args []string) {
		s := cobrashell.New()
		// display welcome info.
		s.Println(fmt.Sprintf("Welcome to Goalign Console %s", version.Version))
		s.Println("type \"help\" to get a list of available commands")
		cobrashell.AddCommands(s, cmd.Root(), nil, cmd.Root().Commands()...)
		// We open a gotree console to interactively execute commands
		s.Run()
	},
}

// Read sequences (possibly not aligned) from a fasta file
func readsequences(file string) (sequences align.SeqBag, err error) {
	var fi goio.Closer
	var r *bufio.Reader

	if fi, r, err = utils.GetReader(file); err != nil {
		return
	}
	defer fi.Close()

	if sequences, err = fasta.NewParser(r).ParseUnalign(); err != nil {
		return
	}

	return
}

// Read aligned sequences from an input file
func readalign(file string) (alchan align.AlignChannel, err error) {
	var fi goio.Closer
	var r *bufio.Reader
	var format int

	if fi, r, err = utils.GetReader(file); err != nil {
		return
	}
	if rootAutoDetectInputFormat {
		if alchan, format, err = utils.ParseMultiAlignmentsAuto(fi, r, rootinputstrict); err != nil {
			return
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
				if err2 := phylip.NewParser(r, rootinputstrict).ParseMultiple(alchan.Achan); err2 != nil {
					alchan.Err = err2
				}
				fi.Close()
			}()
		} else if rootnexus {
			var al align.Alignment
			if al, err = nexus.NewParser(r).Parse(); err != nil {
				return
			}
			alchan.Achan <- al
			fi.Close()
			close(alchan.Achan)
		} else if rootclustal {
			var al align.Alignment
			if al, err = clustal.NewParser(r).Parse(); err != nil {
				return
			}
			alchan.Achan <- al
			fi.Close()
			close(alchan.Achan)
		} else {
			var al align.Alignment
			if al, err = fasta.NewParser(r).Parse(); err != nil {
				return
			}
			alchan.Achan <- al
			fi.Close()
			close(alchan.Achan)
		}
	}
	return
}

// Execute adds all child commands to the root command sets flags appropriately.
// This is called by main.main(). It only needs to happen once to the rootCmd.
func Execute() {
	if err := RootCmd.Execute(); err != nil {
		fmt.Println(err)
		os.Exit(1)
	}
}

func init() {
	cobra.OnInitialize(initConfig)

	RootCmd.PersistentFlags().StringVarP(&infile, "align", "i", "stdin", "Alignment input file")
	RootCmd.PersistentFlags().BoolVarP(&rootphylip, "phylip", "p", false, "Alignment is in phylip? default fasta")
	RootCmd.PersistentFlags().BoolVarP(&rootnexus, "nexus", "x", false, "Alignment is in nexus? default fasta")
	RootCmd.PersistentFlags().BoolVarP(&rootclustal, "clustal", "u", false, "Alignment is in clustal? default fasta")
	RootCmd.PersistentFlags().IntVarP(&rootcpus, "threads", "t", 1, "Number of threads")

	RootCmd.PersistentFlags().Int64Var(&seed, "seed", -1, "Random Seed: -1 = nano seconds since 1970/01/01 00:00:00")
	RootCmd.PersistentFlags().BoolVar(&rootinputstrict, "input-strict", false, "Strict phylip input format (only used with -p)")
	RootCmd.PersistentFlags().BoolVar(&rootoutputstrict, "output-strict", false, "Strict phylip output format (only used with -p)")
	RootCmd.PersistentFlags().BoolVar(&rootoutputoneline, "one-line", false, "Write Phylip sequences on 1 line (only used with -p)")
	RootCmd.PersistentFlags().BoolVar(&rootoutputnoblock, "no-block", false, "Write Phylip sequences without space separated blocks (only used with -p)")

	RootCmd.PersistentFlags().BoolVar(&rootAutoDetectInputFormat, "auto-detect", false, "Auto detects input format (overrides -p, -x and -u)")

	RootCmd.SetHelpTemplate(helptemplate)
}

// initConfig reads in config file and ENV variables if set.
func initConfig() {

}

func writeAlign(al align.Alignment, f *os.File) {
	if rootphylip {
		f.WriteString(phylip.WriteAlignment(al, rootoutputstrict, rootoutputoneline, rootoutputnoblock))
	} else if rootnexus {
		f.WriteString(nexus.WriteAlignment(al))
	} else if rootclustal {
		f.WriteString(clustal.WriteAlignment(al))
	} else {
		f.WriteString(fasta.WriteAlignment(al))
	}
}

func writeAlignString(al align.Alignment) (out string) {
	if rootphylip {
		out = phylip.WriteAlignment(al, rootoutputstrict, rootoutputoneline, rootoutputnoblock)
	} else if rootnexus {
		out = nexus.WriteAlignment(al)
	} else if rootclustal {
		out = clustal.WriteAlignment(al)
	} else {
		out = fasta.WriteAlignment(al)
	}
	return
}

func alignExtension() (out string) {
	if rootphylip {
		out = ".ph"
	} else if rootnexus {
		out = ".nx"
	} else if rootclustal {
		out = ".clustal"
	} else {
		out = ".fa"
	}
	return
}

func writeSequences(seqs align.SeqBag, f *os.File) {
	f.WriteString(fasta.WriteAlignment(seqs))
}

func writeAlignFasta(al align.Alignment, f *os.File) {
	f.WriteString(fasta.WriteAlignment(al))
}

func writeAlignPhylip(al align.Alignment, f *os.File) {
	f.WriteString(phylip.WriteAlignment(al, rootoutputstrict, rootoutputoneline, rootoutputnoblock))
}

func writeAlignNexus(al align.Alignment, f *os.File) {
	f.WriteString(nexus.WriteAlignment(al))
}

func writeAlignClustal(al align.Alignment, f *os.File) {
	f.WriteString(clustal.WriteAlignment(al))
}

func writeAlignPaml(al align.Alignment, f *os.File) {
	f.WriteString(paml.WriteAlignment(al))
}

func openWriteFile(file string) (f *os.File, err error) {
	if file == "stdout" || file == "-" {
		f = os.Stdout
	} else if file == "none" {
		f, err = os.OpenFile(os.DevNull, os.O_APPEND|os.O_CREATE|os.O_WRONLY, 0666)
	} else {
		f, err = os.Create(file)
	}
	return
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

func closeWriteFile(f goio.Closer, filename string) {
	if filename != "-" && filename != "stdout" && filename != "none" {
		f.Close()
	}
}
