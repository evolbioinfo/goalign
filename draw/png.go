package draw

import (
	"bufio"
	"image"
	"image/color"
	"image/png"

	"github.com/evolbioinfo/goalign/align"
)

// Colors from bioSyntax (doi.org/10.1186/s12859-018-2315-y)
var nucleotideColors = map[rune]color.RGBA{
	'A': {71, 255, 25, 255}, 'C': {255, 70, 65, 255}, 'G': {240, 144, 0, 255}, 'T': {65, 146, 255, 255},
	'R': {255, 254, 128, 255}, 'Y': {225, 128, 255, 255}, 'S': {255, 155, 128, 255}, 'W': {128, 255, 242, 255},
	'K': {144, 184, 44, 255}, 'M': {206, 136, 52, 255}, 'B': {248, 193, 192, 255}, 'D': {199, 255, 185, 255},
	'H': {191, 216, 249, 255}, 'V': {255, 227, 185, 255}, 'N': {230, 230, 230, 255}, '-': {255, 255, 255, 255},
	'X': {230, 230, 230, 255}, '.': {255, 255, 255, 255},
}

// Colors adapted from "Shapely" color scheme http://acces.ens-lyon.fr/biotic/rastop/help/colour.htm
var aminoAcidColors = map[rune]color.RGBA{
	'A': {140, 255, 140, 255}, 'G': {255, 255, 255, 255}, 'L': {69, 94, 69, 255}, 'S': {255, 112, 66, 255},
	'V': {255, 140, 255, 255}, 'T': {184, 76, 0, 255}, 'K': {71, 71, 184, 255}, 'D': {160, 0, 66, 255},
	'I': {0, 76, 0, 255}, 'N': {255, 124, 112, 255}, 'E': {102, 0, 0, 255}, 'P': {82, 82, 82, 255},
	'R': {0, 0, 124, 255}, 'F': {83, 76, 66, 255}, 'Q': {255, 76, 76, 255}, 'Y': {140, 112, 76, 255},
	'H': {112, 112, 255, 255}, 'C': {255, 255, 112, 255}, 'M': {184, 160, 66, 255}, 'W': {79, 70, 0, 255},
	'B': {255, 0, 255, 255}, 'Z': {255, 0, 255, 255}, 'X': {184, 184, 184, 255}, '-': {0, 0, 0, 255},
	'.': {0, 0, 0, 255},
}

type pngLayout struct {
	writer *bufio.Writer
}

func NewPngLayout(writer *bufio.Writer) AlignLayout {
	return &pngLayout{writer}
}

func (layout *pngLayout) DrawAlign(a align.Alignment) (err error) {

	colors := nucleotideColors
	if a.Alphabet() == align.AMINOACIDS {
		colors = aminoAcidColors
	}

	height := a.NbSequences()
	width := a.Length()

	img := image.NewRGBA(image.Rect(0, 0, width, height))

	for j, seq := range a.Sequences() {
		for i, char := range seq.Sequence() {
			c, ok := colors[char]
			if !ok {
				// Set unkown character to black
				c = color.RGBA{0, 0, 0, 0}
			}
			img.Set(i, j, c)
		}
	}

	png.Encode(layout.writer, img)

	return nil
}
