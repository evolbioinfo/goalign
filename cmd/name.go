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
	"github.com/fredericlemoine/goalign/io"
	"github.com/spf13/cobra"
	"os"
)

var trimNb int
var trimMapout string
var trimAlignOut string

// nameCmd represents the name command
var nameCmd = &cobra.Command{
	Use:   "name",
	Short: "This command trims names of sequences",
	Long: `This command trims names of sequences.

It trims sequence names to n characters. It will output mapping between 
old names and new names into a map file as well as the new alignment.
If n is > than seq name, it will add 000 at the end.
If n is < than seq name, it will trim. Names may be identical after this step.
In this case, it will add a unique identifier to the identical names.

Example of usage:

goalign trim name -i align.phy -p -n 10 -m map.txt

`,
	Run: func(cmd *cobra.Command, args []string) {
		if namemap, err := rootalign.TrimNames(trimNb); err != nil {
			io.ExitWithMessage(err)
		} else {
			writeAlign(rootalign, trimAlignOut)
			writeNameMap(namemap, trimMapout)
		}
	},
}

func writeNameMap(namemap map[string]string, outfile string) {
	f, err := os.Create(outfile)
	if err != nil {
		io.ExitWithMessage(err)
	}
	for long, short := range namemap {
		f.WriteString(long)
		f.WriteString("\t")
		f.WriteString(short)
		f.WriteString("\n")
	}
	f.Close()
}

func init() {
	trimCmd.AddCommand(nameCmd)

	nameCmd.PersistentFlags().IntVarP(&trimNb, "nb-char", "n", 1, "Number of characters to keep in sequence names")
	nameCmd.PersistentFlags().StringVarP(&trimMapout, "out-map", "m", "none", "Mapping output file")
	nameCmd.PersistentFlags().StringVarP(&trimAlignOut, "out-align", "o", "stdout", "Renamed alignment output file")

}
