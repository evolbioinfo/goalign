package distance

import (
	"github.com/fredericlemoine/goalign/align"
)

type PDistModel struct {
}

func NewPDistModel() *PDistModel {
	return &PDistModel{}
}

/* computes p-distance between 2 sequences */
func (m *PDistModel) Distance(seq1 []rune, seq2 []rune, weights []float64) float64 {
	diff := countDiffs(seq1, seq2, weights)
	diff = diff / float64(len(seq1))
	return diff
}

func (m *PDistModel) InitModel(al align.Alignment) {

}
