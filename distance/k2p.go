package distance

import (
	"github.com/fredericlemoine/goalign/align"
	"math"
)

type K2PModel struct {
}

func NewK2PModel() *K2PModel {
	return &K2PModel{}
}

/* computes K2P distance between 2 sequences */
func (m *K2PModel) Distance(seq1 []rune, seq2 []rune, weights []float64) float64 {
	trS, trV := countMutations(seq1, seq2, weights)
	trS, trV = trS/float64(len(seq1)), trV/float64(len(seq1))
	dist := -0.5*math.Log(1-2*trS-trV) - 0.25*math.Log(1-2*trV)
	if dist > 0 {
		return (dist)
	} else {
		return (0)
	}
}

func (m *K2PModel) InitModel(al align.Alignment) {

}
