package cmd

import (
	"fmt"
	"image/color"

	"github.com/spf13/cobra"
	"gonum.org/v1/plot"
	"gonum.org/v1/plot/font"
	"gonum.org/v1/plot/palette/brewer"
	"gonum.org/v1/plot/plotter"
	"gonum.org/v1/plot/vg"
	"gonum.org/v1/plot/vg/draw"

	"github.com/evolbioinfo/goalign/align"
	"github.com/evolbioinfo/goalign/distance/dna"
	"github.com/evolbioinfo/goalign/io"
	"github.com/evolbioinfo/goalign/io/utils"
)

var simplotOutput string
var simplotModel string
var simplotRefSeq string
var simplotWindows int
var simplotWindowStep int
var simplotOutImageFile string
var simplotOutImageWidth int
var simplotOutImageHeight int
var simplotGroup bool
var simplotGroupSep string
var simplotGroupField int

// computedistCmd represents the computedist command
var computeSimplotCmd = &cobra.Command{
	Use:   "simplot",
	Short: "Compute simplot data",
	Long: `Compute simplot data

`,
	RunE: func(cmd *cobra.Command, args []string) (err error) {
		var f utils.StringWriterCloser
		var aligns *align.AlignChannel
		var p *plot.Plot
		var point *plotter.Scatter
		var line *plotter.Line

		var windows []struct {
			WindowStart, WindowEnd int
			CompSeq                string
			Distance               float64
		}

		if f, err = utils.OpenWriteFile(simplotOutput); err != nil {
			io.LogError(err)
			return
		}
		defer utils.CloseWriteFile(f, computedistOutput)

		if aligns, err = readalign(infile); err != nil {
			io.LogError(err)
			return
		}

		fmt.Fprintf(f, "start\tend\tcomp\tdist\n")
		for align := range aligns.Achan {

			if aligns.Err != nil {
				err = aligns.Err
				io.LogError(err)
			}
			if windows, err = dna.SimPlotDistances(align, simplotRefSeq,
				simplotModel, simplotWindows, simplotWindowStep, simplotGroup,
				simplotGroupSep, simplotGroupField); err != nil {
				io.LogError(err)
				return
			}

			for _, w := range windows {
				fmt.Fprintf(f, "%d\t%d\t%s\t%v\n", w.WindowStart, w.WindowEnd, w.CompSeq, w.Distance)
			}

			if simplotOutImageFile != "none" {
				p = plot.New()
				p.Title.Text = "Simplot"
				p.X.Label.Text = "Genomic position"
				p.Y.Label.Text = "Similarity to reference"

				colors := brewer.DivergingPalette{
					ID:         "div",
					Name:       "div",
					Laptop:     brewer.Good,
					CRT:        brewer.Good,
					ColorBlind: brewer.Good,
					Copy:       brewer.Good,
					Projector:  brewer.Good,
					Color: []color.Color{
						color.RGBA{R: 212, G: 222, B: 69, A: 255},
						color.RGBA{R: 116, G: 209, B: 94, A: 255},
						color.RGBA{R: 67, G: 89, B: 186, A: 255},
						color.RGBA{R: 199, G: 76, B: 184, A: 255},
						color.RGBA{R: 247, G: 184, B: 13, A: 255},
					},
				}

				pts := make(map[string]plotter.XYs)
				for _, l := range windows {
					xy := plotter.XY{X: float64(l.WindowStart), Y: 1.0 - l.Distance}
					pts[l.CompSeq] = append(pts[l.CompSeq], xy)
				}

				i := 0
				t := plotter.PaletteThumbnailers(colors)
				for k := range pts {
					line, point, err = plotter.NewLinePoints(pts[k])
					fmt.Println(line)
					point.Shape = draw.CircleGlyph{}
					point.Radius = 1
					point.Color = colors.Colors()[i%len(colors.Colors())]
					line.Color = colors.Colors()[i%len(colors.Colors())]
					line.StepStyle = plotter.NoStep
					line.Width = 2
					p.Add(line)  //, point)
					p.Add(point) //, point)
					p.Legend.Add(k, t[i%len(colors.Colors())])
					i++
				}
				// Save the plot to a file.
				if err = p.Save(font.Length(simplotOutImageWidth)*vg.Inch, font.Length(simplotOutImageHeight)*vg.Inch, simplotOutImageFile); err != nil {
					io.LogError(err)
					return
				}
			}

		}
		return
	},
}

func init() {
	computeCmd.AddCommand(computeSimplotCmd)
	computeSimplotCmd.PersistentFlags().StringVarP(&simplotOutput, "output", "o", "stdout", "Distance matrix output file")
	computeSimplotCmd.PersistentFlags().StringVarP(&simplotModel, "model", "m", "k2p", "Model for distance computation")
	computeSimplotCmd.PersistentFlags().StringVarP(&simplotRefSeq, "refseq", "r", "-", "Reference sequence to compare all others")
	computeSimplotCmd.PersistentFlags().IntVarP(&simplotWindows, "window-size", "w", 100, "Window size")
	computeSimplotCmd.PersistentFlags().IntVarP(&simplotWindowStep, "window-step", "s", 100, "Window step")
	computeSimplotCmd.PersistentFlags().StringVar(&simplotOutImageFile, "image", "none", "LTT plot image image output file")
	computeSimplotCmd.PersistentFlags().IntVar(&simplotOutImageWidth, "image-width", 4, "LTT plot image image output width")
	computeSimplotCmd.PersistentFlags().IntVar(&simplotOutImageHeight, "image-height", 4, "LTT plot image output heigh")
	computeSimplotCmd.PersistentFlags().BoolVar(&simplotGroup, "group", false, "If sequences must be grouped")
	computeSimplotCmd.PersistentFlags().IntVar(&simplotGroupField, "field", 0, "Field number for extracting group (if group is true)")
	computeSimplotCmd.PersistentFlags().StringVar(&simplotGroupSep, "sep", "_", "Separator for extracting group (if group is true)")
}
