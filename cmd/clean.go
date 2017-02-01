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

var cleanOutput string
var allgaps bool

// cleanCmd represents the clean command
var cleanCmd = &cobra.Command{
	Use:   "clean",
	Short: "Removes gap sites",
	Long: `Removes sites constituted of gaps

If -a is given, then removes sites having only gaps
if -a is not given, removes sites having at least one gap.

`,
	Run: func(cmd *cobra.Command, args []string) {
		f := openWriteFile(cleanOutput)
		for al := range rootaligns {
			al.RemoveGaps(allgaps)
			writeAlign(al, f)
		}
		f.Close()

	},
}

func init() {
	RootCmd.AddCommand(cleanCmd)
	cleanCmd.PersistentFlags().BoolVarP(&allgaps, "all-gaps", "a", false, "Remove positions having only gaps")
	cleanCmd.PersistentFlags().StringVarP(&cleanOutput, "output", "o", "stdout", "Cleaned alignment output file")
}
