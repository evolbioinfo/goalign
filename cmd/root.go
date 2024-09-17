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
	"strconv"
	"strings"
	"time"

	"github.com/evolbioinfo/goalign/align"
	"github.com/evolbioinfo/goalign/io/clustal"
	"github.com/evolbioinfo/goalign/io/fasta"
	"github.com/evolbioinfo/goalign/io/nexus"
	"github.com/evolbioinfo/goalign/io/paml"
	"github.com/evolbioinfo/goalign/io/partition"
	"github.com/evolbioinfo/goalign/io/phylip"
	"github.com/evolbioinfo/goalign/io/stockholm"
	"github.com/evolbioinfo/goalign/io/utils"
	"github.com/evolbioinfo/goalign/version"
	"github.com/fredericlemoine/cobrashell"
	"github.com/spf13/cobra"
)

var infile string
var rootphylip bool
var rootnexus bool
var rootclustal bool
var rootstockholm bool
var rootcpus int
var rootinputstrict bool = false
var rootoutputstrict bool = false
var rootoutputoneline = false
var rootoutputnoblock = false
var rootalphabet string
var rootAutoDetectInputFormat bool
var seed int64 = -1
var unaligned bool
var ignoreidentical = align.IGNORE_NONE

var helptemplate string = `{{with or .Long .Short }}{{. | trim}}

{{end}}Version: ` + version.Version + `

{{if or .Runnable .HasSubCommands}}{{.UsageString}}{{end}}

If you use the Gotree/Goalign toolkit, please cite:
Lemoine F, Gascuel O. 
Gotree/Goalign: toolkit and Go API to facilitate the development of phylogenetic workflows. 
NAR Genom Bioinform. 2021 Aug 11;3(3):lqab075.
doi: 10.1093/nargab/lqab075. PMID: 34396097; PMCID: PMC8356961.
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

	p := fasta.NewParser(r)
	p.IgnoreIdentical(ignoreidentical)
	if sequences, err = p.ParseUnalign(); err != nil {
		return
	}

	return
}

// Read aligned sequences from an input file
func readalign(file string) (alchan *align.AlignChannel, err error) {
	var fi goio.Closer
	var r *bufio.Reader
	var format int
	var alphabet int
	alchan = &align.AlignChannel{}

	switch rootalphabet {
	case "auto":
		alphabet = align.BOTH
	case "nt":
		alphabet = align.NUCLEOTIDS
	case "aa":
		alphabet = align.AMINOACIDS
	default:
		err = fmt.Errorf("given alphabet is not supported: %s", rootalphabet)
		return
	}

	if fi, r, err = utils.GetReader(file); err != nil {
		return
	}
	if rootAutoDetectInputFormat {
		if alchan, format, err = utils.ParseMultiAlignmentsAuto(fi, r, rootinputstrict, alphabet); err != nil {
			return
		}
		if format == align.FORMAT_PHYLIP {
			rootphylip = true
		} else if format == align.FORMAT_NEXUS {
			rootnexus = true
		} else if format == align.FORMAT_CLUSTAL {
			rootclustal = true
		} else if format == align.FORMAT_STOCKHOLM {
			rootstockholm = true
		}
	} else {
		if rootphylip {
			alchan.Achan = make(chan align.Alignment, 15)
			go func() {
				pp := phylip.NewParser(r, rootinputstrict)
				pp.Alphabet(alphabet)
				pp.IgnoreIdentical(ignoreidentical)
				pp.ParseMultiple(alchan)
				fi.Close()
			}()
		} else if rootnexus {
			var al align.Alignment
			np := nexus.NewParser(r)
			np.Alphabet(alphabet)
			np.IgnoreIdentical(ignoreidentical)
			if al, err = np.Parse(); err != nil {
				return
			}
			alchan.Achan = make(chan align.Alignment, 1)
			alchan.Achan <- al
			fi.Close()
			close(alchan.Achan)
		} else if rootclustal {
			var al align.Alignment
			cp := clustal.NewParser(r)
			cp.Alphabet(alphabet)
			cp.IgnoreIdentical(ignoreidentical)
			if al, err = cp.Parse(); err != nil {
				return
			}
			alchan.Achan = make(chan align.Alignment, 1)
			alchan.Achan <- al
			fi.Close()
			close(alchan.Achan)
		} else if rootstockholm {
			var al align.Alignment
			cp := stockholm.NewParser(r)
			cp.Alphabet(alphabet)
			cp.IgnoreIdentical(ignoreidentical)
			if al, err = cp.Parse(); err != nil {
				return
			}
			alchan.Achan = make(chan align.Alignment, 1)
			alchan.Achan <- al
			fi.Close()
			close(alchan.Achan)
		} else {
			var al align.Alignment
			fp := fasta.NewParser(r)
			fp.Alphabet(alphabet)
			fp.IgnoreIdentical(ignoreidentical)
			if al, err = fp.Parse(); err != nil {
				return
			}
			alchan.Achan = make(chan align.Alignment, 1)
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
	RootCmd.PersistentFlags().BoolVarP(&rootstockholm, "stockholm", "k", false, "Alignment is in stockholm? default fasta")
	RootCmd.PersistentFlags().IntVarP(&rootcpus, "threads", "t", 1, "Number of threads")

	// If ignore is IGNORE_NONE: Does not ignore anything
	// If ignore is IGNORE_NAME: Ignore sequences having the same name (keep the first one whatever their sequence)
	// If ignore is IGNORE_SEQUENCE: Ignore sequences having the same name and the same sequence
	// Otherwise, sets IGNORE_NONE

	RootCmd.PersistentFlags().IntVar(&ignoreidentical, "ignore-identical", align.IGNORE_NONE, fmt.Sprintf("Ignore duplicated sequences that have the same name and potentially have same sequences, %d : Does not ignore anything, %d: Ignore sequences having the same name (keep the first one whatever their sequence), %d: Ignore sequences having the same name and the same sequence", align.IGNORE_NONE, align.IGNORE_NAME, align.IGNORE_SEQUENCE))
	RootCmd.PersistentFlags().Int64Var(&seed, "seed", -1, "Random Seed: -1 = nano seconds since 1970/01/01 00:00:00")
	RootCmd.PersistentFlags().BoolVar(&rootinputstrict, "input-strict", false, "Strict phylip input format (only used with -p)")
	RootCmd.PersistentFlags().BoolVar(&rootoutputstrict, "output-strict", false, "Strict phylip output format (only used with -p)")
	RootCmd.PersistentFlags().BoolVar(&rootoutputoneline, "one-line", false, "Write Phylip sequences on 1 line (only used with -p)")
	RootCmd.PersistentFlags().BoolVar(&rootoutputnoblock, "no-block", false, "Write Phylip sequences without space separated blocks (only used with -p)")
	RootCmd.PersistentFlags().StringVar(&rootalphabet, "alphabet", "auto", "Alignment/Sequences alphabet: auto (default), aa, or nt")

	RootCmd.PersistentFlags().BoolVar(&rootAutoDetectInputFormat, "auto-detect", false, "Auto detects input format (overrides -p, -x and -u)")

	RootCmd.SetHelpTemplate(helptemplate)
}

// initConfig reads in config file and ENV variables if set.
func initConfig() {

}

func writeAlign(al align.Alignment, f utils.StringWriterCloser) {
	if rootphylip {
		f.WriteString(phylip.WriteAlignment(al, rootoutputstrict, rootoutputoneline, rootoutputnoblock))
	} else if rootnexus {
		f.WriteString(nexus.WriteAlignment(al))
	} else if rootclustal {
		f.WriteString(clustal.WriteAlignment(al))
	} else if rootstockholm {
		f.WriteString(stockholm.WriteAlignment(al))
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
	} else if rootstockholm {
		out = ".sto"
	} else {
		out = ".fa"
	}
	return
}

func writeSequences(seqs align.SeqBag, f utils.StringWriterCloser) {
	f.WriteString(fasta.WriteAlignment(seqs))
}

func writeAlignFasta(al align.Alignment, f utils.StringWriterCloser) {
	f.WriteString(fasta.WriteAlignment(al))
}

func writeAlignPhylip(al align.Alignment, f utils.StringWriterCloser) {
	f.WriteString(phylip.WriteAlignment(al, rootoutputstrict, rootoutputoneline, rootoutputnoblock))
}

func writeAlignNexus(al align.Alignment, f utils.StringWriterCloser) {
	f.WriteString(nexus.WriteAlignment(al))
}

func writeAlignClustal(al align.Alignment, f utils.StringWriterCloser) {
	f.WriteString(clustal.WriteAlignment(al))
}

func writeAlignPaml(al align.Alignment, f utils.StringWriterCloser) {
	f.WriteString(paml.WriteAlignment(al))
}

func readMapFile(file string, revert bool) (map[string]string, error) {
	outmap := make(map[string]string)
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
	line, e := utils.Readln(reader)
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
		line, e = utils.Readln(reader)
		nl++
	}

	if err = mapfile.Close(); err != nil {
		return outmap, err
	}

	return outmap, nil
}

func parsePartition(partitionfile string, alilength int) (ps *align.PartitionSet, err error) {
	var f goio.Closer
	var r *bufio.Reader

	if f, r, err = utils.GetReader(partitionfile); err != nil {
		return
	}
	defer f.Close()
	p := partition.NewParser(r)
	ps, err = p.Parse(alilength)
	return
}

func parseIntFile(file string) (ints []int, err error) {
	var f *os.File
	var r *bufio.Reader
	var gr *gzip.Reader
	var c int

	ints = make([]int, 0)

	if file == "stdin" || file == "-" {
		f = os.Stdin
	} else {
		if f, err = os.Open(file); err != nil {
			return
		}
	}

	if strings.HasSuffix(file, ".gz") {
		if gr, err = gzip.NewReader(f); err != nil {
			return
		}
		r = bufio.NewReader(gr)
	} else {
		r = bufio.NewReader(f)
	}
	l, e := utils.Readln(r)
	for e == nil {
		if c, err = strconv.Atoi(l); err != nil {
			return
		}
		ints = append(ints, c)
		l, e = utils.Readln(r)
	}

	return
}
