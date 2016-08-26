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
	"github.com/spf13/cobra"
)

// trimCmd represents the trim command
var trimCmd = &cobra.Command{
	Use:   "trim",
	Short: "This command trims names of sequences or sequences themselves",
	Long: `This command trims names of sequences or sequences themselves.

With "names" subcommand, you can trim names to n characters. In this case, it will
also output mapping between old names and new names into a map file as well as the
new alignment.

With "seq" subcommand, you can trim sequences from start or from end, by n characters.
`,
}

func init() {
	RootCmd.AddCommand(trimCmd)
}
