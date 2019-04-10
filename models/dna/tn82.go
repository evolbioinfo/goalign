package dna

import (
	"github.com/evolbioinfo/goalign/align"
	"math"
)

type TN82Model struct {
	/* Vector of nt proba */
	pi            []float64 // proba of each nt
	numSites      float64   // Number of selected sites (no gaps)
	selectedSites []bool    // true for selected sites
	removegaps    bool      // If true, we will remove posision with >=1 gaps
}

func NewTN82Model(removegaps bool) *TN82Model {
	return &TN82Model{
		nil,
		0,
		nil,
		removegaps,
	}
}

/* computes TN82 distance between 2 sequences */
func (m *TN82Model) Distance(seq1 []rune, seq2 []rune, weights []float64) (float64, error) {
	diff, total := countDiffs(seq1, seq2, m.selectedSites, weights)
	diff = diff / total

	psi := init2DFloat(4, 4)
	totalPairs, err := countNtPairs2Seq(seq1, seq2, m.selectedSites, weights, psi)
	if err != nil {
		return 0.0, err
	}
	for i := 0; i < 4; i++ {
		for j := 0; j < 4; j++ {
			psi[i][j] = psi[i][j] / totalPairs
		}
	}
	denom := 0.0
	for i := 0; i < 4; i++ {
		for j := i + 1; j < 4; j++ {
			denom += psi[i][j] * psi[i][j] / (2 * m.pi[i] * m.pi[j])
		}
	}
	b1 := diff * diff / denom
	dist := -1.0 * b1 * math.Log(1.0-diff/b1)
	if dist > 0 {
		return dist, nil
	} else {
		return 0, nil
	}
}

func (m *TN82Model) InitModel(al align.Alignment, weights []float64) (err error) {
	m.numSites, m.selectedSites = selectedSites(al, weights, m.removegaps)
	m.pi, err = probaNt(al, m.selectedSites, weights)
	return
}
