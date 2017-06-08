// +build ignore

package cmd

import (
	"fmt"
	"github.com/fredericlemoine/goalign/distance"
	"github.com/spf13/cobra"
	"math/rand"
	"os"
	"time"
)

// tntweightbootCmd represents the tntweightboot command
var tntweightbootCmd = &cobra.Command{
	Use:   "tntweightboot",
	Short: "Generates a data file to give as TNT input, with continuous weight for sites",
	Long: `Generates a data file to give as TNT input, with continuous weight for sites

Weights follow a Dirichlet distribution D(n;1,...,1)

`,
	Run: func(cmd *cobra.Command, args []string) {
		rand.Seed(weightbootSeed)

		al := <-rootaligns

		mintnt := 1.0
		maxtnt := 1000.0
		for i := 0; i < weightbootnb; i++ {
			var f *os.File
			if weightbootOutput != "stdout" && weightbootOutput != "-" {
				f = openWriteFile(fmt.Sprintf("%s%d.tnt", weightbootOutput, i))
			} else {
				f = openWriteFile(weightbootOutput)
			}
			var weights []float64 = nil
			weights = distance.BuildWeightsDirichlet(al)
			minw := -1.0
			maxw := -1.0
			for _, w := range weights {
				if w > maxw || maxw == -1 {
					maxw = w
				}
				if w < minw || minw == -1 {
					minw = w
				}
			}
			f.WriteString("xread\n\n")
			f.WriteString("'Tnt input file'\n\n")
			f.WriteString(fmt.Sprintf("%d %d\n", al.Length(), al.NbSequences()))
			al.Iterate(func(name string, sequence string) {
				f.WriteString(fmt.Sprintf("%s %s\n", name, sequence))
			})
			f.WriteString("ccode\n")
			for i, w := range weights {
				tntweight := (w-minw)/(maxw-minw)*(maxtnt-mintnt) + (mintnt)
				if i > 0 {
					f.WriteString(" ")
				}
				f.WriteString(fmt.Sprintf("/%d-[%d", int(tntweight), i))
			}
			f.WriteString("\n")
			f.WriteString(";\n")
			f.Close()
		}
	},
}

func init() {
	buildCmd.AddCommand(tntweightbootCmd)
	tntweightbootCmd.PersistentFlags().Int64VarP(&weightbootSeed, "seed", "s", time.Now().UTC().UnixNano(), "Initial Random Seed")
	tntweightbootCmd.PersistentFlags().StringVarP(&weightbootOutput, "output-prefix", "o", "stdout", "Prefix of Tnt input data file")
	tntweightbootCmd.PersistentFlags().IntVarP(&weightbootnb, "nboot", "n", 1, "Number of bootstrap replicates to build")
}
