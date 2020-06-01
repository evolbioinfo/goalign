package cmd

import (
	"fmt"
	"os"
	"path/filepath"
	"strconv"

	"github.com/evolbioinfo/goalign/align"
	"github.com/evolbioinfo/goalign/io"
	"github.com/spf13/cobra"
)

var subsitesout string = "stdout"
var subsitesfile string
var subsitesrefseq string
var subsitesreverse bool

// subsitesCmd represents the subsites command
var subsitesCmd = &cobra.Command{
	Use:   "subsites",
	Short: "Takes a subset of the sites of the alignment",
	Long: `Takes a subset of the sites of the alignment

It takes an alignment and extracts only some sites from it, given
a set of sites given on the command line or on a given file.
If some sites are outside the alignment, will exit with an error.

For example:
goalign subsites -p -i al.phy 1 2 3 4

will select sites 1, 2, 3, 4 from the alignment (0-based inclusive)

The output format is the same than input format.

If --ref-seq <name> is specified, then the coordinates are considered according to the 
given sequence, and without considering gaps.

For example:
If al.fa is:
>s1
--ACG--AT-GC
>s2
GGACGTTATCGC

goalign subsites -i al.fa --ref-seq s1 1 2 3

will output:

>s1
CGA
>s2
CGA

If --reverse is given, then will output all positions but the ones that should be output.
`,
	RunE: func(cmd *cobra.Command, args []string) (err error) {
		var aligns *align.AlignChannel
		var f *os.File
		var c int
		var subalign align.Alignment

		if aligns, err = readalign(infile); err != nil {
			io.LogError(err)
			return
		}
		if f, err = openWriteFile(subsitesout); err != nil {
			io.LogError(err)
			return
		}

		refseq := cmd.Flags().Changed("ref-seq")

		var sites []int

		if subsitesfile != "none" {
			if sites, err = parseIntFile(subsitesfile); err != nil {
				io.LogError(err)
				return
			}
		} else {
			sites = make([]int, 0)
			for _, s := range args {
				if c, err = strconv.Atoi(s); err != nil {
					io.LogError(err)
					return
				}
				sites = append(sites, c)
			}
		}

		if len(sites) == 0 {
			err = fmt.Errorf("No sites are provided")
			io.LogError(err)
			return
		}

		fileid := ""
		filenum := 0
		extension := filepath.Ext(subsitesout)
		name := subsitesout[0 : len(subsitesout)-len(extension)]

		for al := range aligns.Achan {

			if filenum > 0 && subsitesout != "stdout" && subsitesout != "-" {
				fileid = fmt.Sprintf("%s_al%d%s", name, filenum, extension)
				f.Close()
				if f, err = openWriteFile(fileid); err != nil {
					io.LogError(err)
					return
				}
			}

			sitespos := sites
			if refseq {
				sitespos, err = al.RefSites(subsitesrefseq, sites)
			}

			if subsitesreverse {
				if sitespos, err = al.InversePositions(sitespos); err != nil {
					io.LogError(err)
					return
				}
			}
			if subalign, err = al.SelectSites(sitespos); err != nil {
				io.LogError(err)
				return
			}
			writeAlign(subalign, f)
			filenum++
		}

		f.Close()

		if aligns.Err != nil {
			err = aligns.Err
			io.LogError(err)
		}
		return
	},
}

func init() {
	RootCmd.AddCommand(subsitesCmd)
	subsitesCmd.PersistentFlags().StringVarP(&subsitesout, "output", "o", "stdout", "Alignment output file")
	subsitesCmd.PersistentFlags().StringVar(&subsitesrefseq, "ref-seq", "none", "Reference sequence on which coordinates are given")
	subsitesCmd.PersistentFlags().StringVar(&subsitesfile, "sitefile", "none", "File with positions of sites to select (one perline)")
	subsitesCmd.PersistentFlags().BoolVarP(&subsitesreverse, "reverse", "r", false, "Take all but the given sites")
}
