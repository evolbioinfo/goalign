package cmd

import (
	"bufio"
	"fmt"
	"path/filepath"

	"github.com/evolbioinfo/goalign/align"
	"github.com/evolbioinfo/goalign/draw"
	"github.com/evolbioinfo/goalign/io"
	"github.com/evolbioinfo/goalign/io/utils"
	"github.com/spf13/cobra"
)

// pngCmd represents the png command
var pngCmd = &cobra.Command{
	Use:   "png",
	Short: "Draw alignments in a png file",
	Long: `Draw alignments in a png file

One line per sequence, one pixel per character. 
Color schemes are specific to the alphabet of sequences: 
	- The nucleotide colors are from bioSyntax (doi.org/10.1186/s12859-018-2315-y).
	- The amino acid colors are adapted from "Shapely colours" 
	  (http://acces.ens-lyon.fr/biotic/rastop/help/colour.htm)
`,
	RunE: func(cmd *cobra.Command, args []string) (err error) {
		var l draw.AlignLayout
		var aligns *align.AlignChannel
		var f utils.StringWriterCloser

		if aligns, err = readalign(infile); err != nil {
			io.LogError(err)
			return
		}

		nalign := 0
		for al := range aligns.Achan {
			fname := drawOutput
			// Add an index to file output name
			// if there are several alignments to draw
			if nalign > 0 {
				ext := filepath.Ext(fname)
				fname = fmt.Sprintf("%s_%d.%s", fname[0:len(fname)-len(ext)], nalign, ext)
			}
			if f, err = utils.OpenWriteFile(fname); err != nil {
				io.LogError(err)
				return
			}
			al.CleanNames(nil)
			w := bufio.NewWriter(f)
			l = draw.NewPngLayout(w)
			l.DrawAlign(al)
			w.Flush()
			f.Close()
			nalign++
		}

		if aligns.Err != nil {
			err = aligns.Err
			io.LogError(err)
		}
		return
	},
}

func init() {
	drawCmd.AddCommand(pngCmd)
}
