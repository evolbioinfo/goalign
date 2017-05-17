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
	"fmt"

	"github.com/fredericlemoine/goalign/io"
	"github.com/spf13/cobra"
)

var cleanOutput string
var cleanCutoff float64
var cleanQuiet bool

// cleanCmd represents the clean command
var cleanCmd = &cobra.Command{
	Use:   "clean",
	Short: "Removes gap sites",
	Long: `Removes sites constituted of gaps

Removes alignments sites constitued of >= cutoff gap sites.

Exception for a cutoff of 0: removes sites constitued of > 0 gap sites.

Examples:
- With a cutoff of 0.5: a site with 5 gaps over 10 sequences will be removed;
- With a cutoff of 0.5: a site with 4 gaps over 10 sequences will not be removed;
- With a cutoff of 0.0 a site with 1 gap over 10 sequences will be removed.

If cutoff is <0 or >1, it will be considered as 0, which means that every site with at least 1 gap
will be removed.
`,
	Run: func(cmd *cobra.Command, args []string) {
		f := openWriteFile(cleanOutput)
		i := 0
		for al := range rootaligns {
			beforelength := al.Length()
			al.RemoveGaps(cleanCutoff)
			afterlength := al.Length()
			writeAlign(al, f)
			if !cleanQuiet {
				io.PrintMessage(fmt.Sprintf("Alignment (%d) length before cleaning=%d", i, beforelength))
				io.PrintMessage(fmt.Sprintf("Alignment (%d) length after cleaning=%d", i, afterlength))
				io.PrintMessage(fmt.Sprintf("Alignment (%d) number of gaps=%d", i, beforelength-afterlength))
			}
		}
		f.Close()

	},
}

func init() {
	RootCmd.AddCommand(cleanCmd)
	cleanCmd.PersistentFlags().StringVarP(&cleanOutput, "output", "o", "stdout", "Cleaned alignment output file")
	cleanCmd.PersistentFlags().Float64VarP(&cleanCutoff, "cutoff", "c", 0, "Cutoff for gap deletion : 0 remove sites with > 0 gap, 1 remove sites with 100% gaps)")
	cleanCmd.PersistentFlags().BoolVarP(&cleanQuiet, "quiet", "q", false, "Do not print results on stderr")
}
