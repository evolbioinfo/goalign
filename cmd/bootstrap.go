package cmd

import (
	"archive/tar"
	"bufio"
	"compress/gzip"
	"errors"
	"fmt"
	"os"
	"time"

	"github.com/evolbioinfo/goalign/align"
	"github.com/evolbioinfo/goalign/io"
	"github.com/spf13/cobra"
)

var bootstrapNb int
var bootstrapoutprefix string
var bootstrapOrder bool
var bootstraptar bool
var bootstrapgz bool
var bootstrapfrac float64
var bootstrappartitionstr string
var bootstrapoutputpartitionstr string

// seqbootCmd represents the bootstrap command
var seqbootCmd = &cobra.Command{
	Use:   "seqboot",
	Short: "Generates n bootstrap alignments from an input alignment",
	Long: `Generates n bootstrap alignments from input alignment.

The input may be a Phylip or Fasta file.

- With -S option, sequence order is shuffled in output alignments.
- With --tar-gz option, only one .tar.gz file is generated, containing 
  all alignments files.

- If the input file is in phylip format (-p option), then output bootstrap
  files will be in phylip format as well.
- If the input file is in fasta format (no -p option), then output bootstrap
  files will be in fasta format as well.

- It is possible to give a initial seed (--seed). In this case several runs of 
  the tool will give the exact same results.

- If frac is < 1.0, output alignments are partial bootstrap alignments as is phylip
  seqboot, which means that the sites are sampled from the full alignment with 
  replacement, but the bootstrap alignment length is a fraction of the original alignment.

Example of usage:

goalign build seqboot -i align.phylip -p -n 500 -o boot --tar-gz
goalign build seqboot -i align.phylip -p -n 500 -o boot_ 
`,
	RunE: func(cmd *cobra.Command, args []string) (err error) {
		var alignChan *align.AlignChannel
		var aligns []align.Alignment
		var al align.Alignment
		var f *os.File
		var tw *tar.Writer
		var gw *gzip.Writer
		var inputpartition, outputpartition *align.PartitionSet
		var bootstring string
		var boot, tmpboot align.Alignment

		// We read input alignment
		if alignChan, err = readalign(infile); err != nil {
			io.LogError(err)
			return
		}

		if bootstrapoutprefix == "none" {
			err = errors.New("Output bootstrap file prefix is mandatory")
			io.LogError(err)
			return
		}

		// We take the first alignment of the channel
		al, _ = <-alignChan.Achan
		if alignChan.Err != nil {
			err = alignChan.Err
			io.LogError(err)
			return
		}

		// If a partition file is given, then we parse it
		if bootstrappartitionstr != "none" {
			if inputpartition, err = parsePartition(bootstrappartitionstr, al.Length()); err != nil {
				io.LogError(err)
				return
			}
			if err = inputpartition.CheckSites(); err != nil {
				io.LogError(err)
				return
			}
			if aligns, err = al.Split(inputpartition); err != nil {
				io.LogError(err)
				return
			}
			outputpartition = align.NewPartitionSet(al.Length())
			//fmt.Println(bootstrappartition.String())
		} else {
			aligns = []align.Alignment{al}
		}

		// Create new tar(/gz) file
		if bootstraptar {
			if bootstrapgz {
				f, err = os.Create(bootstrapoutprefix + ".tar.gz")
			} else {
				f, err = os.Create(bootstrapoutprefix + ".tar")
			}
			if err != nil {
				io.LogError(err)
				return
			}
			defer f.Close()
			if bootstrapgz {
				gw = gzip.NewWriter(f)
				defer gw.Close()
				tw = tar.NewWriter(gw)
			} else {
				tw = tar.NewWriter(f)
			}
			defer tw.Close()
		}

		for idx := 0; idx < bootstrapNb; idx++ {
			boot = nil
			bootid := bootstrapoutprefix + fmt.Sprintf("%d", idx)
			// There may be several alignments to process if there are
			// several partitions. We generate bootstrap replicates
			// for each partition, and then concatenate them all.
			for _, a := range aligns {
				tmpboot = a.BuildBootstrap(bootstrapfrac)
				if boot == nil {
					boot = tmpboot
				} else {
					if err = boot.Concat(tmpboot); err != nil {
						io.LogError(err)
						return
					}
				}
			}
			// We shuffle sequence order
			if bootstrapOrder {
				boot.ShuffleSequences()
			}

			bootstring = writeAlignString(boot)

			// Output
			if bootstraptar {
				if err = addstringtotargz(tw, bootid, bootstring); err != nil {
					io.LogError(err)
					return
				}
			} else {
				if err = writenewfile(bootid+alignExtension(), bootstrapgz, bootstring); err != nil {
					io.LogError(err)
					return
				}
			}
		}

		var start, end int = 0, 0
		if outputpartition != nil {
			for i, a := range aligns {
				start = end
				end = start + a.Length()
				// We initialize an outputpartition
				// Which will have all the sites of each
				// partition grouped together.
				outputpartition.AddRange(
					inputpartition.PartitionName(i),
					inputpartition.ModeleName(i),
					start, end-1, 1)
			}
			if bootstrapoutputpartitionstr == "" {
				bootstrapoutputpartitionstr = bootstrappartitionstr + "_boot"
			}
			writenewfile(bootstrapoutputpartitionstr, false, outputpartition.String())
		}

		return
	},
}

func writenewfile(name string, gz bool, bootstring string) (err error) {
	var f *os.File

	if gz {
		if f, err = os.Create(name + ".gz"); err != nil {
			return
		}
		gw := gzip.NewWriter(f)
		buf := bufio.NewWriter(gw)
		buf.WriteString(bootstring)
		buf.Flush()
		gw.Close()
		f.Close()
	} else {
		if f, err = os.Create(name); err != nil {
			return
		}
		f.WriteString(bootstring)
		f.Close()
	}
	return
}

func addstringtotargz(tw *tar.Writer, name string, align string) error {
	// now lets create the header as needed for this file within the tarball
	alignbytes := []byte(align)
	ext := alignExtension()
	header := new(tar.Header)
	header.Name = name + ext
	header.Size = int64(len(alignbytes))
	header.Mode = 436 //int64(stat.Mode())
	header.ModTime = time.Now()
	// write the header to the tarball archive
	if err := tw.WriteHeader(header); err != nil {
		return err
	}

	// copy the file data to the tarball
	if _, err := tw.Write(alignbytes); err != nil {
		return err
	}
	return nil
}

func minInt(a, b int) int {
	if a < b {
		return a
	}
	return b
}

func init() {
	buildCmd.AddCommand(seqbootCmd)

	seqbootCmd.PersistentFlags().BoolVarP(&bootstrapOrder, "shuf-order", "S", false, "Also shuffle order of sequences in bootstrap files")
	seqbootCmd.PersistentFlags().BoolVar(&bootstraptar, "tar", false, "Will create a single tar file with all bootstrap alignments (one thread for tar, but not a bottleneck)")
	seqbootCmd.PersistentFlags().BoolVar(&bootstrapgz, "gz", false, "Will gzip output file(s). Maybe slow if combined with --tar (only one thread working for tar/gz)")
	seqbootCmd.PersistentFlags().IntVarP(&bootstrapNb, "nboot", "n", 1, "Number of bootstrap replicates to build")
	seqbootCmd.PersistentFlags().Float64VarP(&bootstrapfrac, "frac", "f", 1.0, "Fraction of sites to sample (if < 1.0: Partial bootstrap as in phylip seqboot)")
	seqbootCmd.PersistentFlags().StringVar(&bootstrappartitionstr, "partition", "none", "File containing definition of the partitions")
	seqbootCmd.PersistentFlags().StringVar(&bootstrapoutputpartitionstr, "out-partition", "", "File containing output partitions (default: same name as input partition with _boot suffix)")
	seqbootCmd.PersistentFlags().StringVarP(&bootstrapoutprefix, "out-prefix", "o", "none", "Prefix of output bootstrap files")
}
