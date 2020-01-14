package cmd

import (
	"fmt"
	"os"

	"github.com/spf13/cobra"

	"github.com/evolbioinfo/goalign/align"
	"github.com/evolbioinfo/goalign/io"
)

var splitpartition *align.PartitionSet
var splitpartitionstr string
var splitprefix string

// seqbootCmd represents the bootstrap command
var splitCmd = &cobra.Command{
	Use:   "split",
	Short: "Splits an input alignment given a partition file",
	Long: `Splits an input alignment given a partition file.

Output alignment files will be in the same format as input alignment, 
with file names corresponding to partition names.

Example of usage:
goalign split -i align.phylip --partition partition.txt 
`,
	RunE: func(cmd *cobra.Command, args []string) (err error) {
		var aligns *align.AlignChannel
		var splitAligns []align.Alignment

		var f *os.File

		if aligns, err = readalign(infile); err != nil {
			io.LogError(err)
			return
		}

		align, _ := <-aligns.Achan
		if aligns.Err != nil {
			err = aligns.Err
			io.LogError(err)
			return
		}

		if splitpartitionstr != "none" {
			if splitpartition, err = parsePartition(splitpartitionstr, align.Length()); err != nil {
				io.LogError(err)
				return
			}
			if err = splitpartition.CheckSites(); err != nil {
				io.LogError(err)
				return
			}
		} else {
			err = fmt.Errorf("Partition file must be provided")
			io.LogError(err)
			return
		}

		if splitAligns, err = align.Split(splitpartition); err != nil {
			io.LogError(err)
			return
		}

		for i, a := range splitAligns {
			name := splitprefix + splitpartition.PartitionName(i) + alignExtension()
			if f, err = openWriteFile(name); err != nil {
				io.LogError(err)
				return
			}
			writeAlign(a, f)
			f.Close()
		}

		return
	},
}

func init() {
	RootCmd.AddCommand(splitCmd)

	splitCmd.PersistentFlags().StringVarP(&splitprefix, "out-prefix", "o", "", "Prefix of output files")
	splitCmd.PersistentFlags().StringVar(&splitpartitionstr, "partition", "none", "File containing definition of the partitions")
}
