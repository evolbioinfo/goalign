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
	"github.com/spf13/cobra"
)

var rogueNb float64
var rogueNameFile string

// rogueCmd represents the rogue command
var rogueCmd = &cobra.Command{
	Use:   "rogue",
	Short: "Simulate rogue taxa",
	Long:  `Simulate rogue `,
	Run: func(cmd *cobra.Command, args []string) {
		f := openWriteFile(shuffleOutput)
		nameFile := openWriteFile(rogueNameFile)
		for al := range rootaligns {
			names, _ := al.SimulateRogue(rogueNb)
			writeAlign(al, f)
			for _, n := range names {
				nameFile.WriteString(n)
				nameFile.WriteString("\n")
			}
		}
		nameFile.Close()
		f.Close()
	},
}

func init() {
	shuffleCmd.AddCommand(rogueCmd)

	rogueCmd.PersistentFlags().Float64VarP(&rogueNb, "prop-seq", "n", 0.5, "Proportion of the  sequences to consider as rogue")
	rogueCmd.PersistentFlags().StringVar(&rogueNameFile, "rogue-file", "stdout", "Rogue sequence names output file")
}
