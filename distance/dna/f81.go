package dna

import (
	"math"

	"github.com/evolbioinfo/goalign/align"
)

type F81Model struct {
	pi            []float64 // Vector of nt stationary proba
	b1            float64   // Parameter for distance computation
	numSites      float64   // Number of selected sites (no gaps)
	selectedSites []bool    // true for selected sites
	removegaps    bool      // If true, we will remove posision with >=1 gaps
	gamma         bool
	alpha         float64
}

func NewF81Model(removegaps bool) *F81Model {
	return &F81Model{
		nil,
		0,
		0,
		nil,
		removegaps,
		false,
		0.,
	}
}

/* computes F81 distance between 2 sequences */
func (m *F81Model) Distance(seq1 []rune, seq2 []rune, weights []float64) (float64, error) {
	var dist float64

	diff, total := countDiffs(seq1, seq2, m.selectedSites, weights)
	diff = diff / total

	if m.gamma {
		dist = 1. * m.b1 * m.alpha * (math.Pow(1.-diff/m.b1, -1./m.alpha) - 1.)
	} else {
		dist = -1. * m.b1 * math.Log(1.-diff/m.b1)
	}
	if dist > 0 {
		return dist, nil
	} else {
		return 0, nil
	}
}

func (m *F81Model) InitModel(al align.Alignment, weights []float64, gamma bool, alpha float64) (err error) {
	m.gamma = gamma
	m.alpha = alpha
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
