package cmd

import (
	"fmt"

	"github.com/evolbioinfo/goalign/align"
	"github.com/evolbioinfo/goalign/io"

	"github.com/evolbioinfo/goalign/io/utils"
	"github.com/spf13/cobra"
)

var convertgffrefseq string    // reference sequence on which the coordinates are given
var convertgffdestseq string   // sequence on which the coordinates should be converted
var convertgffcoordfile string // file with all coordinates of the sequences to extract
var convertgffoutput string

// extractCmd represents the extract command
var convertGFFCmd = &cobra.Command{
	Use:   "convertgff",
	Short: "Converts the coordinates of an input GFF file from one reference sequence to another",
	Long: `This command converts the coordinates of an input GFF file from one reference sequence to another.
The input GFF file should be specified with --coordinates, and the reference sequence on which the coordinates are given should be specified with --ref-seq.
The output file should be specified with --output.
Example:
goalign convertgff --gff input.gff --ref-seq ref1 --dest-seq dest1 --output output.gff
`,
	RunE: func(cmd *cobra.Command, args []string) (err error) {
		var aligns *align.AlignChannel
		var f utils.StringWriterCloser
		var subcoords []extractSubSequence
		var newstart, newlength int

		if convertgffcoordfile == "none" {
			err = fmt.Errorf("subsequence gff file should be specified")
			return
		}

		if aligns, err = readalign(infile); err != nil {
			io.LogError(err)
			return
		}

		if subcoords, err = parseGFFFile(convertgffcoordfile); err != nil {
			io.LogError(err)
			return
		}

		al := <-aligns.Achan
		if aligns.Err != nil {
			err = aligns.Err
			io.LogError(err)
			return
		}

		if f, err = utils.OpenWriteFile(convertgffoutput); err != nil {
			io.LogError(err)
			return
		}

		for _, subseq := range subcoords {
			//fmt.Printf("Log: Converting coordinates for %s: [%v,%v[ on %s\n", subseq.name, subseq.starts, subseq.ends, subseq.chromosome)
			for i, s := range subseq.starts {
				e := subseq.ends[i]
				l := e - s
				//fmt.Printf("Log: Converting coordinates for %s: [%d,%d[ on %s\n", subseq.name, s, e, subseq.chromosome)
				if s < 0 || e > al.Length() {
					err = fmt.Errorf("coordinates are outside alignment: [%d,%d[", s, e)
					io.LogError(err)
					return
				}
				if s >= e {
					err = fmt.Errorf("block length should be >0 : [%d,%d[", s, e)
					io.LogError(err)
					return
				}

				if s, l, err = al.RefCoordinates(convertgffrefseq, s, l); err != nil {
					io.LogError(err)
					return
				}
				if newstart, newlength, err = al.DestCoordinates(convertgffdestseq, s, l); err != nil {
					io.LogError(err)
					return
				}
				fmt.Fprintf(f, "%s\t%s\t%s\t%d\t%d\t.\t%s\t.\tID=gene-%s;Name=%s;gbkey=%s;gene=%s;gene_biotype=%s\n", convertgffdestseq, subseq.source, "gene", newstart+1, newstart+newlength, strand2string(subseq.strand), subseq.name, subseq.name, "Gene", subseq.name, "protein_coding")
				fmt.Fprintf(f, "%s\t%s\t%s\t%d\t%d\t.\t%s\t.\tID=cds-%s;Name=%s;gbkey=%s;Parent=gene-%s\n", convertgffdestseq, subseq.source, "CDS", newstart+1, newstart+newlength, strand2string(subseq.strand), subseq.name, subseq.name, "CDS", subseq.name)

			}
		}
		f.Close()

		return
	},
}

func init() {
	RootCmd.AddCommand(convertGFFCmd)
	convertGFFCmd.PersistentFlags().StringVar(&convertgffrefseq, "ref-seq", "none", "Reference sequence on which coordinates are given")
	convertGFFCmd.PersistentFlags().StringVarP(&convertgffoutput, "output", "o", "stdout", "Output folder")
	convertGFFCmd.PersistentFlags().StringVar(&convertgffcoordfile, "gff", "none", "GFF file with all coordinates to convert")
	convertGFFCmd.PersistentFlags().StringVar(&convertgffdestseq, "dest-seq", "none", "Destination sequence to which coordinates are converted")

}
