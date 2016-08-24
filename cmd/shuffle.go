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
	"github.com/fredericlemoine/goalign/io/fasta"
	"github.com/fredericlemoine/goalign/io/phylip"
	"github.com/spf13/cobra"
	"math/rand"
	"os"
	"time"
)

var shuffleSeed int64
var shuffleOutput string

// shuffleCmd represents the shuffle command
var shuffleCmd = &cobra.Command{
	Use:   "shuffle",
	Short: "A brief description of your command",
	Long: `A longer description that spans multiple lines and likely contains examples
and usage of using your command. For example:

Cobra is a CLI library for Go that empowers applications.
This application is a tool to generate the needed files
to quickly create a Cobra application.`,
	PersistentPreRun: func(cmd *cobra.Command, args []string) {
		RootCmd.PersistentPreRun(cmd, args)
		rand.Seed(shuffleSeed)
	},
	PersistentPostRun: func(cmd *cobra.Command, args []string) {
		var f *os.File
		var err error

		if shuffleOutput == "stdout" {
			f = os.Stdout
		} else {
			f, err = os.Create(shuffleOutput)
			if err != nil {
				panic(err)
			}
		}
		if rootphylip {
			f.WriteString(phylip.WriteAlignment(rootalign))
		} else {
			f.WriteString(fasta.WriteAlignment(rootalign))
		}

		f.Close()
	},
}

func init() {
	RootCmd.AddCommand(shuffleCmd)

	shuffleCmd.PersistentFlags().Int64VarP(&shuffleSeed, "seed", "s", time.Now().UTC().UnixNano(), "Initial Random Seed")
	shuffleCmd.PersistentFlags().StringVarP(&shuffleOutput, "output", "o", "stdout", "Shuffled alignment output file")
}
