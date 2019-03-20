package cmd

import (
	"fmt"
	"os"
	"path/filepath"

	"github.com/fredericlemoine/goalign/align"
	"github.com/fredericlemoine/goalign/io"
	"github.com/spf13/cobra"
)

var subseqout string = "stdout"
var subseqstart int
var subseqlength int
var subseqstep int

// subseqCmd represents the subseq command
var subseqCmd = &cobra.Command{
	Use:   "subseq",
	Short: "Take a sub-alignment from the input alignment",
	Long: `Take a sub-alignment from the input alignment

It takes an alignment and extracts sub-sequences from it, given
a start position (0-based inclusive) and a length.
If the length is after the end of the alignment, will stop at the 
end of the alignment.

For example:
goalign subseq -p -i al.phy -s 9 -l 10

This will extract a sub-alignment going from 10th position, with a length of 10.

The output format is the same than input format.

Sliding window:
---------------
If --step is given and > 0, then Several sub-alignments will be produced,
and corresponding to all alignments in windows of sizes -l, and with starts:
[start, start+step, ..., end-length].

Example with an alignment al.phy of size 10 (0123456789)

goalign subseq -i al.phy -s 0 -l 5 --step 1 will produce:

01234
 12345
  23456
   34567
    45678
     56789

Warning: If output is stdout, it works only if input format is Phylip, because 
it is possible to split alignments afterwards (goalign divide for example).

Several alignments:
------------------
If several alignments are present in the input file and the output is a file 
(not stdout or -) , then :
* First alignment, first subalignment: results will be placed in the given file
  (ex out.fasta)
* First alignment, other subalignments (sliding windows): results will be placed
  in file with the given name with "_sub<i>" suffix (ex: out_sub1.fasta, out_sub2.fasta, etc.)
* Other alignments, first subalignment: results will be placed in the given file
  with "_al<i>" suffix (ex out_al1.fasta, out_al2.fasta, etc.)
* Other alignments, other subalignments: results will be placed in the given file
  with "_al<i>" and "_sub<i> suffixes (ex out_al1_sub1.fasta, out_al1_sub2.fasta, etc.)
`,
	RunE: func(cmd *cobra.Command, args []string) (err error) {
		var aligns *align.AlignChannel
		var f *os.File
		var subalign align.Alignment

		if aligns, err = readalign(infile); err != nil {
			io.LogError(err)
			return
		}
		if f, err = openWriteFile(subseqout); err != nil {
			io.LogError(err)
			return
		}

		fileid := ""
		filenum := 0
		extension := filepath.Ext(subseqout)
		name := subseqout[0 : len(subseqout)-len(extension)]

		for al := range aligns.Achan {
			start := subseqstart
			if filenum > 0 && subseqout != "stdout" && subseqout != "-" {
				fileid = fmt.Sprintf("_al%d", filenum)
				f.Close()
				if f, err = openWriteFile(name + fileid + extension); err != nil {
					io.LogError(err)
					return
				}
			}
			subalignnum := 0
			for {
				if subalign, err = al.SubAlign(start, subseqlength); err != nil {
					io.LogError(err)
					return
				}
				writeAlign(subalign, f)
				start += subseqstep
				if subseqstep == 0 || (start+subseqlength) > al.Length() {
					break
				} else {
					if subseqout != "stdout" && subseqout != "-" {
						subalignnum++
						f.Close()
						n := fmt.Sprintf("%s%s_sub%d%s", name, fileid, subalignnum, extension)
						if f, err = openWriteFile(n); err != nil {
							io.LogError(err)
							return
						}
					}
				}
			}
			filenum++
		}
		f.Close()

		if aligns.Err != nil {
			err = aligns.Err
			io.LogError(err)
		}
		return
	},
}

func init() {
	RootCmd.AddCommand(subseqCmd)
	subseqCmd.PersistentFlags().StringVarP(&subseqout, "output", "o", "stdout", "Alignment output file")
	subseqCmd.PersistentFlags().IntVarP(&subseqstart, "start", "s", 0, "Start position (0-based inclusive)")
	subseqCmd.PersistentFlags().IntVarP(&subseqlength, "length", "l", 10, "Length of the sub alignment")
	subseqCmd.PersistentFlags().IntVar(&subseqstep, "step", 0, "Step: If > 0, then will generate several alignments, for each window of length l, with starts: [start,start+step, ..., end-l]* ")
}
