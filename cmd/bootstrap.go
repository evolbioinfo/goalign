package cmd

import (
	"archive/tar"
	"compress/gzip"
	"fmt"
	"github.com/fredericlemoine/goalign/io/fasta"
	"github.com/fredericlemoine/goalign/io/phylip"
	"github.com/spf13/cobra"
	"math/rand"
	"os"
	"time"
)

var bootstrapSeed int64
var bootstrapNb int
var bootstrapoutprefix string
var bootstrapOrder bool
var bootstraptargz bool

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
		var f *os.File
		var err error
		var tw *tar.Writer
		var gw *gzip.Writer
		var bootstring string

		if bootstraptargz {
			f, err = os.Create(bootstrapoutprefix + ".tar.gz")
			if err != nil {
				panic(err)
			}
			defer f.Close()
			gw = gzip.NewWriter(f)
			defer gw.Close()
			tw = tar.NewWriter(gw)
			defer tw.Close()
		}

		for i := 0; i < bootstrapNb; i++ {
			bootid := bootstrapoutprefix + fmt.Sprintf("%d", i)
			boot := rootalign.BuildBootstrap()
			if bootstrapOrder {
				boot.ShuffleSequences()
			}

			if rootphylip {
				bootstring = phylip.WriteAlignment(boot)
			} else {
				bootstring = fasta.WriteAlignment(boot)
			}

			if bootstraptargz {
				if err = addstringtotargz(tw, bootid, bootstring); err != nil {
					panic(err)
				}
			} else {
				f, err = os.Create(bootid)
				if err != nil {
					panic(err)
				}
				f.WriteString(bootstring)
				f.Close()
			}
		}
	},
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
	bootstrapCmd.PersistentFlags().BoolVar(&bootstraptargz, "tar-gz", false, "Will create a single targz file with all bootstrap alignments")
	bootstrapCmd.PersistentFlags().IntVarP(&bootstrapNb, "nboot", "n", 1, "Number of bootstrap replicates to build")
	bootstrapCmd.PersistentFlags().StringVarP(&bootstrapoutprefix, "out-prefix", "o", "none", "Prefix of output bootstrap files")
}
