// Copyright © 2016 NAME HERE <EMAIL ADDRESS>
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

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
	Short: "A brief description of your command",
	Long: `A longer description that spans multiple lines and likely contains examples
and usage of using your command. For example:

Cobra is a CLI library for Go that empowers applications.
This application is a tool to generate the needed files
to quickly create a Cobra application.`,
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
	bootstrapCmd.PersistentFlags().IntVarP(&bootstrapNb, "nboot", "n", 1, "Prefix of output bootstrap files")
	bootstrapCmd.PersistentFlags().StringVarP(&bootstrapoutprefix, "prefix", "r", "none", "Number of bootstrap replicates to build")
}
