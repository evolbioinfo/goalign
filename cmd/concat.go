package cmd

import (
	"fmt"
	"path/filepath"
	"strings"

	"github.com/evolbioinfo/goalign/align"
	"github.com/evolbioinfo/goalign/io"
	"github.com/evolbioinfo/goalign/io/utils"
	"github.com/spf13/cobra"
)

var concatout string
var concatlog string
var concatpart string

// concatCmd represents the concat command
var concatCmd = &cobra.Command{
	Use:   "concat",
	Short: "Concatenates a set of alignments",
	Long: `Concatenates a set of alignments.
For example:

If format is phylip, it may contain several alignments in one file. 
Then we can concatenate all of them:
goalign concat -i align.phy

If format is Fasta, it is not possible, then you must give other alignments in the form:
goalign concat -i align.fasta others*.fasta

It is possible to give only otherfiles, without -i, by giving -i none
   goalign concat -i none align*.fasta
or goalign concat -i none -p align*.phy

If -l is given, the coordinates of all the input alignments in the concatenated alignment 
are written in the log file (tab separated values : start (0-based inclusive), end (0-based exclusive),
input file name).

If --out-partition is provided, a partition file is written at the specified location, 
with each partition matching its source alignment. The default model is GTR for nucleotide 
alignments and LG for aminoacid alignments.

`,
	RunE: func(cmd *cobra.Command, args []string) (err error) {
		var aligns *align.AlignChannel
		var alchan *align.AlignChannel
		var outputpartition *align.PartitionSet
		var partitionmodel string = "GTR"
		var curalign align.Alignment
		var f utils.StringWriterCloser // output align
		var p utils.StringWriterCloser // output partition
		var l utils.StringWriterCloser // output log
		var start int
		if l, err = utils.OpenWriteFile(concatlog); err != nil {
			io.LogError(err)
			return
		}

		start = 0
		if infile != "none" {
			if aligns, err = readalign(infile); err != nil {
				io.LogError(err)
				return
			}
			name := strings.TrimSuffix(filepath.Base(infile), filepath.Ext(infile))
			nal := 0
			for al := range aligns.Achan {
				fmt.Fprintf(l, "%d\t%d\t%s\n", start, start+al.Length(), infile)
				if curalign == nil {
					curalign = al
					outputpartition = align.NewPartitionSet(al.Length())
					if al.Alphabet() == align.AMINOACIDS {
						partitionmodel = "LG"
					}
				} else {
					if err = curalign.Concat(al); err != nil {
						io.LogError(err)
						return
					}
					outputpartition.ExtendAliLength(al.Length())
				}
				partname := name
				if nal > 0 {
					partname = fmt.Sprintf("%s%d", name, nal)
				}
				if err = outputpartition.AddRange(partname, "Model", start, start+al.Length()-1, 1); err != nil {
					io.LogError(err)
					return
				}
				start += al.Length()
				nal++
			}
			if aligns.Err != nil {
				err = aligns.Err
				io.LogError(err)
				return
			}
		}

		for _, otherfile := range args {
			if alchan, err = readalign(otherfile); err != nil {
				io.LogError(err)
				return
			}
			name := strings.TrimSuffix(filepath.Base(otherfile), filepath.Ext(otherfile))
			nal := 0
			for al := range alchan.Achan {
				fmt.Fprintf(l, "%d\t%d\t%s\n", start, start+al.Length(), otherfile)
				if curalign == nil {
					curalign = al
					outputpartition = align.NewPartitionSet(al.Length())
					if al.Alphabet() == align.AMINOACIDS {
						partitionmodel = "LG"
					}
				} else {
					if err = curalign.Concat(al); err != nil {
						io.LogError(err)
						return
					}
					outputpartition.ExtendAliLength(al.Length())
				}
				partname := name
				if nal > 0 {
					partname = fmt.Sprintf("%s%d", name, nal)
				}
				if err = outputpartition.AddRange(partname, partitionmodel, start, start+al.Length()-1, 1); err != nil {
					io.LogError(err)
					return
				}
				start += al.Length()
				nal++
			}
			if alchan.Err != nil {
				err = alchan.Err
				io.LogError(err)
				return
			}
		}
		utils.CloseWriteFile(l, concatlog)

		if err = outputpartition.CheckSites(); err != nil {
			io.LogError(err)
			return
		}
		if p, err = utils.OpenWriteFile(concatpart); err != nil {
			io.LogError(err)
			return
		}
		fmt.Fprintf(p, "%s", outputpartition.String())
		utils.CloseWriteFile(p, concatpart)

		if f, err = utils.OpenWriteFile(concatout); err != nil {
			io.LogError(err)
			return
		}
		writeAlign(curalign, f)
		utils.CloseWriteFile(f, concatout)

		return
	},
}

func init() {
	RootCmd.AddCommand(concatCmd)
	concatCmd.PersistentFlags().StringVarP(&concatout, "output", "o", "stdout", "Alignment output file")
	concatCmd.PersistentFlags().StringVarP(&concatlog, "log", "l", "none", "Log output file (coordinates of all input alignments in the concatenated alignment)")
	concatCmd.PersistentFlags().StringVar(&concatpart, "out-partition", "none", "File containing output partitions")
}
