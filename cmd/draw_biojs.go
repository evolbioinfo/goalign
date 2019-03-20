package cmd

import (
	"bufio"
	"fmt"
	"os"
	"path/filepath"

	"github.com/fredericlemoine/goalign/align"
	"github.com/fredericlemoine/goalign/draw"
	"github.com/fredericlemoine/goalign/io"
	"github.com/spf13/cobra"
)

// pngCmd represents the png command
var biojsCmd = &cobra.Command{
	Use:   "biojs",
	Short: "Draw alignments in html file using msaviewer from biojs",
	Long: `Draw alignments in html file using msaviewer from biojs

See http://msa.biojs.net/ for more informations
`,
	RunE: func(cmd *cobra.Command, args []string) (err error) {
		var l draw.AlignLayout
		var aligns *align.AlignChannel
		var f *os.File

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
			if f, err = openWriteFile(fname); err != nil {
				io.LogError(err)
				return
			}
			al.CleanNames(nil)
			w := bufio.NewWriter(f)
			l = draw.NewBioJSLayout(w)
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
	drawCmd.AddCommand(biojsCmd)
}
