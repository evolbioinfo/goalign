package cmd

import (
	"archive/tar"
	"bufio"
	"compress/gzip"
	"errors"
	"fmt"
	"github.com/spf13/cobra"
	"os"
	"sync"
	"time"

	"github.com/fredericlemoine/goalign/align"
	"github.com/fredericlemoine/goalign/io"
)

var bootstrapNb int
var bootstrapoutprefix string
var bootstrapOrder bool
var bootstraptar bool
var bootstrapgz bool

type outboot struct {
	bootstr string
	name    string
}

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
  the tool will give the exact same results (if run on 1 thread, an thus 
  computations are done on a single thread in this case).

Example of usage:

goalign build seqboot -i align.phylip -p -n 500 -o boot --tar-gz
goalign build seqboot -i align.phylip -p -n 500 -o boot_ 
`,
	RunE: func(cmd *cobra.Command, args []string) (err error) {
		var aligns align.AlignChannel

		if aligns, err = readalign(infile); err != nil {
			io.LogError(err)
			return
		}

		if bootstrapoutprefix == "none" {
			err = errors.New("Output bootstrap file prefix is mandatory")
			io.LogError(err)
			return
		}

		var f *os.File
		var tw *tar.Writer
		var gw *gzip.Writer

		align, _ := <-aligns.Achan
		if aligns.Err != nil {
			err = aligns.Err
			io.LogError(err)
			return
		}

		bootidx := make(chan int, 100)
		outchan := make(chan outboot, 100)

		cpus := rootcpus
		if bootstraptar {
			cpus = min_int(1, cpus-1)
		}

		go func() {
			for i := 0; i < bootstrapNb; i++ {
				bootidx <- i
			}
			close(bootidx)
		}()

		var wg sync.WaitGroup // For waiting end of step computation
		// Seed is set => 1 thread
		if cmd.Flags().Changed("seed") {
			cpus = 1
		}
		for i := 0; i < cpus; i++ {
			wg.Add(1)
			go func() {
				var bootstring string
				for idx := range bootidx {
					bootid := bootstrapoutprefix + fmt.Sprintf("%d", idx)
					boot := align.BuildBootstrap()
					if bootstrapOrder {
						boot.ShuffleSequences()
					}

					bootstring = writeAlignString(boot)

					// Output
					if bootstraptar {
						outchan <- outboot{bootstring, bootid}
					} else {
						if err2 := writenewfile(bootid, bootstrapgz, bootstring); err2 != nil {
							io.LogError(err2)
							err = err2
							return
						}
					}
				}
				wg.Done()
			}()
		}

		go func() {
			wg.Wait()
			close(outchan)
		}()

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

		idx := 0
		for oboot := range outchan {
			if bootstraptar {
				if err = addstringtotargz(tw, oboot.name, oboot.bootstr); err != nil {
					io.LogError(err)
					return
				}
			}
			idx++
		}
		return
	},
}

func writenewfile(name string, gz bool, bootstring string) (err error) {
	var f *os.File

	ext := alignExtension()

	if gz {
		if f, err = os.Create(name + ext + ".gz"); err != nil {
			return
		} else {
			gw := gzip.NewWriter(f)
			buf := bufio.NewWriter(gw)
			buf.WriteString(bootstring)
			buf.Flush()
			gw.Close()
			f.Close()
		}
	} else {
		if f, err = os.Create(name + ext); err != nil {
			return
		} else {
			f.WriteString(bootstring)
			f.Close()
		}
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

func min_int(a, b int) int {
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
	seqbootCmd.PersistentFlags().StringVarP(&bootstrapoutprefix, "out-prefix", "o", "none", "Prefix of output bootstrap files")
}
