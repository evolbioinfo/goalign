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
	"github.com/fredericlemoine/goalign/distance"
	"github.com/spf13/cobra"
)

var k2pgamma bool

// k2pCmd represents the k2p command
var k2pCmd = &cobra.Command{
	Use:   "k2p",
	Short: "Compute k2p distance matrix of 2 sequences",
	Long: `Compute k2p distance matrix of 2 sequences

For example:

goalign computedist k2p -i align.ph -p
goalign computedist k2p -i align.fa

`,
	Run: func(cmd *cobra.Command, args []string) {
		var weights []float64 = nil
		if k2pgamma {
			weights = distance.BuildWeights(rootalign)
		}
		distMatrix = distance.MatrixK2P(rootalign, weights)
	},
}

func init() {
	computedistCmd.AddCommand(k2pCmd)
	k2pCmd.PersistentFlags().BoolVarP(&k2pgamma, "with-gamma", "g", false, "If weights are attributed according to a gamma distribution")

}
