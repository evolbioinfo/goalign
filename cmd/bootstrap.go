package cmd

import (
	"archive/tar"
	"bufio"
	"compress/gzip"
	"fmt"
	"github.com/fredericlemoine/goalign/io/fasta"
	"github.com/fredericlemoine/goalign/io/phylip"
	"github.com/spf13/cobra"
	"math/rand"
	"os"
	"sync"
	"time"
)

var bootstrapSeed int64
var bootstrapNb int
var bootstrapoutprefix string
var bootstrapOrder bool
var bootstraptar bool
var bootstrapgz bool
var bootstrapCpus int

type outboot struct {
	bootstr string
	name    string
}

// bootstrapCmd represents the bootstrap command
var bootstrapCmd = &cobra.Command{
	Use:   "bootstrap",
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

- It is possible to give a initial seed (-s). In this case several runs of 
  the tool will give the exact same results. 

Example of usage:

goalign build bootstrap -i align.phylip -p -n 500 -o boot --tar-gz
goalign build bootstrap -i align.phylip -p -n 500 -o boot_ 
`,
	Run: func(cmd *cobra.Command, args []string) {
		if bootstrapoutprefix == "none" {
			panic("Output bootstrap file prefix is mandatory")
		}

		rand.Seed(bootstrapSeed)
		var err error
		var f *os.File
		var tw *tar.Writer
		var gw *gzip.Writer

		bootidx := make(chan int, 100)
		outchan := make(chan outboot, 100)

		go func() {
			for i := 0; i < bootstrapNb; i++ {
				bootidx <- i
			}
			close(bootidx)
		}()

		var wg sync.WaitGroup // For waiting end of step computation
		for i := 0; i < bootstrapCpus; i++ {
			wg.Add(1)
			go func() {
				var bootstring string
				for idx := range bootidx {
					bootid := bootstrapoutprefix + fmt.Sprintf("%d", idx)
					boot := rootalign.BuildBootstrap()
					if bootstrapOrder {
						boot.ShuffleSequences()
					}
					if rootphylip {
						bootstring = phylip.WriteAlignment(boot)
					} else {
						bootstring = fasta.WriteAlignment(boot)
					}

					// Output
					if bootstraptar {
						outchan <- outboot{bootstring, bootid}
					} else {
						writenewfile(bootid, bootstrapgz, bootstring)
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
				panic(err)
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
					panic(err)
				}
			}
			idx++
		}
	},
}

func writenewfile(name string, gz bool, bootstring string) {
	if gz {
		if f, err := os.Create(name + ".gz"); err != nil {
			panic(err)
		} else {
			gw := gzip.NewWriter(f)
			buf := bufio.NewWriter(gw)
			buf.WriteString(bootstring)
			gw.Close()
			f.Close()
		}
	} else {
		if f, err := os.Create(name); err != nil {
			panic(err)
		} else {
			f.WriteString(bootstring)
			f.Close()
		}
	}
}

func addstringtotargz(tw *tar.Writer, name string, align string) error {
	// now lets create the header as needed for this file within the tarball
	alignbytes := []byte(align)

	header := new(tar.Header)
	header.Name = name
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

func init() {
	buildCmd.AddCommand(bootstrapCmd)

	bootstrapCmd.PersistentFlags().Int64VarP(&bootstrapSeed, "seed", "s", time.Now().UTC().UnixNano(), "Initial Random Seed")
	bootstrapCmd.PersistentFlags().BoolVarP(&bootstrapOrder, "shuf-order", "S", false, "Also shuffle order of sequences in bootstrap files")
	bootstrapCmd.PersistentFlags().BoolVar(&bootstraptar, "tar", false, "Will create a single tar file with all bootstrap alignments (one thread for tar, but not a bottleneck)")
	bootstrapCmd.PersistentFlags().BoolVar(&bootstrapgz, "gz", false, "Will gzip output file(s). Maybe slow if combined with --tar (only one thread working for tar/gz)")
	bootstrapCmd.PersistentFlags().IntVarP(&bootstrapNb, "nboot", "n", 1, "Number of bootstrap replicates to build")
	bootstrapCmd.PersistentFlags().StringVarP(&bootstrapoutprefix, "out-prefix", "o", "none", "Prefix of output bootstrap files")
	bootstrapCmd.PersistentFlags().IntVarP(&bootstrapCpus, "threads", "t", 1, "Number of threads")
}
