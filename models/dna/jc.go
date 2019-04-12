package dna

import (
	"math"

	"github.com/evolbioinfo/goalign/align"
)

type JCModel struct {
	numSites      float64 // Number of selected sites (no gaps)
	selectedSites []bool  // true for selected sites
	removegaps    bool    // If true, we will remove posision with >=1 gaps
}

func NewJCModel(removegaps bool) *JCModel {
	return &JCModel{
		0,
		nil,
		removegaps,
	}
}

/* computes JC69 distance between 2 sequences */
func (m *JCModel) Distance(seq1 []rune, seq2 []rune, weights []float64) (float64, error) {
	diff, total := countDiffs(seq1, seq2, m.selectedSites, weights)
	diff = diff / total
	dist := -3.0 / 4.0 * math.Log(1.0-4.0/3.0*diff)
	if dist > 0 {
		return dist, nil
	} else {
		return 0, nil
	}
}

func (m *JCModel) InitModel(al align.Alignment, weights []float64) (err error) {
	m.numSites, m.selectedSites = selectedSites(al, weights, m.removegaps)
	return
}

// Left vectors and right vectors are given in column-major format
func (m *JCModel) Eigens() (val []float64, leftvectors, rightvectors []float64, err error) {
	val = []float64{
		0,
		-4. / 3.,
		-4. / 3.,
		-4. / 3.,
	}

	leftvectors = []float64{
		1. / 4., 1. / 4., 1. / 4., 1. / 4.,
		-1. / 4., -1. / 4., 3. / 4., -1. / 4.,
		-1. / 4., 3. / 4., -1. / 4., -1. / 4.,
		3. / 4., -1. / 4., -1. / 4., -1. / 4.,
	}

	rightvectors = []float64{
		1., 0., 0., 1.,
		1., 0., 1., 0.,
		1., 1., 0., 0.,
		1., -1., -1., -1.,
	}
	return
}
