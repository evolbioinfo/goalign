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
		for i := 0; i < bootstrapNb; i++ {
			f, err = os.Create(bootstrapoutprefix + fmt.Sprintf("%d", i))
			if err != nil {
				panic(err)
			}

			boot := rootalign.BuildBootstrap()
			if bootstrapOrder {
				boot.ShuffleSequences()
			}
			if rootphylip {
				f.WriteString(phylip.WriteAlignment(boot))
			} else {
				f.WriteString(fasta.WriteAlignment(boot))
			}
			f.Close()
		}
	},
}

func init() {
	buildCmd.AddCommand(bootstrapCmd)

	bootstrapCmd.PersistentFlags().Int64VarP(&bootstrapSeed, "seed", "s", time.Now().UTC().UnixNano(), "Initial Random Seed")
	bootstrapCmd.PersistentFlags().BoolVarP(&bootstrapOrder, "shuf-order", "S", false, "Also shuffle order of sequences in bootstrap files")
	bootstrapCmd.PersistentFlags().IntVarP(&bootstrapNb, "nboot", "n", 1, "Prefix of output bootstrap files")
	bootstrapCmd.PersistentFlags().StringVarP(&bootstrapoutprefix, "prefix", "r", "none", "Number of bootstrap replicates to build")
}
