// Copyright © 2017 NAME HERE <EMAIL ADDRESS>
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
	"os"

	"github.com/fredericlemoine/goalign/align"
	"github.com/fredericlemoine/goalign/io"
	"github.com/spf13/cobra"
)

var rogueNb float64
var rogueLength float64
var rogueNameFile string

// rogueCmd represents the rogue command
var rogueCmd = &cobra.Command{
	Use:   "rogue",
	Short: "Simulate rogue taxa",
	Long: `Simulate rogue by shuffling sites of some sequences.

To do so, it takes a proportion of the sequences and shuffles a proportion l of its sites

Example: We want to simulate 0.25 "rogue taxa" by shuffling 0.25 of their length:

goalign shuffle rogue -i al -n 0.25 -l 0.25

Before        After
S1 12345678    S1 1234567
S2 12345678 => S2 1234567
S3 12345678    S3 1634527 <= This one: 2 nucleotides are shuffled
S4 12345678    S4 1234567
`,
	RunE: func(cmd *cobra.Command, args []string) (err error) {
		var aligns *align.AlignChannel
		var namefile *os.File
		var f *os.File

		if aligns, err = readalign(infile); err != nil {
			io.LogError(err)
			return
		}
		if namefile, err = openWriteFile(rogueNameFile); err != nil {
			io.LogError(err)
			return
		}
		defer closeWriteFile(namefile, rogueNameFile)

		if f, err = openWriteFile(shuffleOutput); err != nil {
			io.LogError(err)
			return
		}
		defer closeWriteFile(f, shuffleOutput)

		for al := range aligns.Achan {
			names, _ := al.SimulateRogue(rogueNb, rogueLength)
			writeAlign(al, f)
			for _, n := range names {
				namefile.WriteString(n)
				namefile.WriteString("\n")
			}
		}

		if aligns.Err != nil {
			err = aligns.Err
			io.LogError(err)
		}
		return
	},
}

func init() {
	shuffleCmd.AddCommand(rogueCmd)

	rogueCmd.PersistentFlags().Float64VarP(&rogueNb, "prop-seq", "n", 0.5, "Proportion of the  sequences to consider as rogue")
	rogueCmd.PersistentFlags().Float64VarP(&rogueLength, "length", "l", 0.5, "Proportion of the sites to shuffle")
	rogueCmd.PersistentFlags().StringVar(&rogueNameFile, "rogue-file", "stdout", "Rogue sequence names output file")
}
