package distance

import (
	"github.com/evolbioinfo/goalign/align"
	"math"
)

type F81Model struct {
	pi            []float64 // Vector of nt stationary proba
	b1            float64   // Parameter for distance computation
	numSites      float64   // Number of selected sites (no gaps)
	selectedSites []bool    // true for selected sites
	removegaps    bool      // If true, we will remove posision with >=1 gaps
}

func NewF81Model(removegaps bool) *F81Model {
	return &F81Model{
		nil,
		0,
		0,
		nil,
		removegaps,
	}
}

/* computes F81 distance between 2 sequences */
func (m *F81Model) Distance(seq1 []rune, seq2 []rune, weights []float64) (float64, error) {
	diff, total := countDiffs(seq1, seq2, m.selectedSites, weights)
	diff = diff / total

	dist := -1.0 * m.b1 * math.Log(1.0-diff/m.b1)
	if dist > 0 {
		return dist, nil
	} else {
		return 0, nil
	}
}

func (m *F81Model) InitModel(al align.Alignment, weights []float64) (err error) {
	m.numSites, m.selectedSites = selectedSites(al, weights, m.removegaps)
	m.b1 = 0.0
	m.pi, err = probaNt(al, m.selectedSites, weights)
	if err == nil {
		for i, _ := range m.pi {
			m.b1 += m.pi[i] * m.pi[i]
		}
		m.b1 = 1 - m.b1
	}
	return
}
