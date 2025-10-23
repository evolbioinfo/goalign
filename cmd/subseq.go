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

var subseqout string = "stdout"
var subseqstart int
var subseqlength int
var subseqstep int
var subseqrefseq string
var subseqreverse bool

// subseqCmd represents the subseq command
var subseqCmd = &cobra.Command{
	Use:   "subseq",
	Short: "Take a sub-alignment from the input alignment",
	Long: `Take a sub-alignment from the input alignment

It takes an alignment and extracts sub-sequences from it, given
a start position (0-based inclusive) and a length.
If the length (l)  is after the end of the alignment, will stop at the 
end of the alignment.
If the length l <0 , then the extracted sequences will be [start,alilength-l[
If the length l <0 and a reference sequence is given, the sub alignment will span [start,reflength-l[
of the ref sequence

For example:
goalign subseq -p -i al.phy -s 9 -l 10

This will extract a sub-alignment going from 10th position, with a length of 10.

The output format is the same than input format.

If --ref-seq <name> is specified, then the coordinates are considered according the 
given sequence, and without considering gaps.

For example:
If al.fa is:
>s1
--ACG--AT-GC
>s2
GGACGTTATCGC

goalign subseq -i al.fa -s 0 -l 4 --ref-seq s1

will output:

>s1
ACG--A
>s2
ACGTTA

Not compatible with --step

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

Not compatible with --ref-seq

Several alignments:
------------------
If several alignments are present in the input file and the output is a file 
(not stdout nor -) , then :
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
		var f utils.StringWriterCloser
		var subalign, subaligntmp align.Alignment

		if aligns, err = readalign(infile); err != nil {
			io.LogError(err)
			return
		}
		if f, err = utils.OpenWriteFile(subseqout); err != nil {
			io.LogError(err)
			return
		}

		fileid := ""
		filenum := 0
		extension := filepath.Ext(subseqout)
		name := subseqout[0 : len(subseqout)-len(extension)]

		refseq := cmd.Flags().Changed("ref-seq")
		if refseq && subseqstep > 0 {
			err = fmt.Errorf("--ref-seq and --step options are not compatible")
			return
		}

		for al := range aligns.Achan {
			start := subseqstart
			leng := subseqlength

			if filenum > 0 && subseqout != "stdout" && subseqout != "-" {
				fileid = fmt.Sprintf("_al%d", filenum)
				f.Close()
				if f, err = utils.OpenWriteFile(name + fileid + extension); err != nil {
					io.LogError(err)
					return
				}
			}
			if refseq {
				if leng < 0 {
					if ref, found := al.GetSequenceChar(subseqrefseq); !found {
						err = fmt.Errorf("reference sequence %s not found in the alignment", subseqrefseq)
						io.LogError(err)
						return
					} else {
						leng = len(strings.ReplaceAll(string(ref), "-", "")) + leng - start
					}
				}
				if start, leng, err = al.RefCoordinates(subseqrefseq, start, leng); err != nil {
					io.LogError(err)
					return
				}
			} else {
				if leng < 0 {
					leng = (al.Length() + leng) - start
				}
			}

			subalignnum := 0
			for {
				starts := []int{start}
				lens := []int{leng}
				if subseqreverse {
					if starts, lens, err = al.InverseCoordinates(start, leng); err != nil {
						io.LogError(err)
						return
					}
				}
				subalign = nil
				for i, s := range starts {
					l := lens[i]
					if subaligntmp, err = al.SubAlign(s, l); err != nil {
						io.LogError(err)
						return
					}
					if subalign == nil {
						subalign = subaligntmp
					} else {
						if err = subalign.Concat(subaligntmp); err != nil {
							io.LogError(err)
							return
						}
					}
				}
				writeAlign(subalign, f)
				start += subseqstep
				if subseqstep == 0 || (start+leng) > al.Length() {
					break
				} else {
					if subseqout != "stdout" && subseqout != "-" {
						subalignnum++
						f.Close()
						n := fmt.Sprintf("%s%s_sub%d%s", name, fileid, subalignnum, extension)
						if f, err = utils.OpenWriteFile(n); err != nil {
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
	subseqCmd.PersistentFlags().IntVarP(&subseqlength, "length", "l", 10, "Length of the sub alignment. If l <0, then the extracted sequences will be [start,alilength-l[")
	subseqCmd.PersistentFlags().StringVar(&subseqrefseq, "ref-seq", "none", "Reference sequence on which coordinates are given")
	subseqCmd.PersistentFlags().BoolVarP(&subseqreverse, "reverse", "r", false, "Take all but the given subsequence")
	subseqCmd.PersistentFlags().IntVar(&subseqstep, "step", 0, "Step: If > 0, then will generate several alignments, for each window of length l, with starts: [start,start+step, ..., end-l]* ")
}
