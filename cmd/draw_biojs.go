package cmd

import (
	"bufio"
	"fmt"
	"path/filepath"

	"github.com/fredericlemoine/goalign/draw"
	"github.com/spf13/cobra"
)

// pngCmd represents the png command
var biojsCmd = &cobra.Command{
	Use:   "biojs",
	Short: "Draw alignments in html file using msaviewer from biojs",
	Long: `Draw alignments in html file using msaviewer from biojs

See http://msa.biojs.net/ for more informations
`,
	Run: func(cmd *cobra.Command, args []string) {
		var l draw.AlignLayout

		nalign := 0
		for al := range rootaligns.Achan {
			fname := drawOutput
			// Add an index to file output name
			// if there are several alignments to draw
			if nalign > 0 {
				ext := filepath.Ext(fname)
				fname = fmt.Sprintf("%s_%d.%s", fname[0:len(fname)-len(ext)], nalign, ext)
			}
			f := openWriteFile(fname)
			w := bufio.NewWriter(f)
			l = draw.NewBioJSLayout(w)
			l.DrawAlign(al)
			w.Flush()
			f.Close()
			nalign++
		}
	},
}

func init() {
	drawCmd.AddCommand(biojsCmd)
}
