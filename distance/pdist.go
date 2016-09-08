package distance

import (
	"github.com/fredericlemoine/goalign/align"
)

type PDistModel struct {
	numSites      float64 // Number of selected sites (no gaps)
	selectedSites []bool  // true for selected sites
}

func NewPDistModel() *PDistModel {
	return &PDistModel{}
}

/* computes p-distance between 2 sequences */
func (m *PDistModel) Distance(seq1 []rune, seq2 []rune, weights []float64) float64 {
	diff, total := countDiffs(seq1, seq2, weights)
	diff = diff / total
	return diff
}

func (m *PDistModel) InitModel(al align.Alignment, weights []float64) {
	m.numSites, m.selectedSites = selectedSites(al, weights)
}
