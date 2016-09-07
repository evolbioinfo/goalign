package distance

import (
	"github.com/fredericlemoine/goalign/align"
	"math"
)

type JCModel struct {
}

func NewJCModel() *JCModel {
	return &JCModel{}
}

/* computes K2P distance between 2 sequences */
func (m *JCModel) Distance(seq1 []rune, seq2 []rune, weights []float64) float64 {
	diff, total := countDiffs(seq1, seq2, weights)
	diff = diff / total
	dist := -3.0 / 4.0 * math.Log(1.0-4.0/3.0*diff)
	if dist > 0 {
		return (dist)
	} else {
		return (0)
	}
}

func (m *JCModel) InitModel(al align.Alignment) {

}
