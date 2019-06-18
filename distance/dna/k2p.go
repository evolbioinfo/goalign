package dna

import (
	"math"

	"github.com/evolbioinfo/goalign/align"
)

type K2PModel struct {
	numSites      float64 // Number of selected sites (no gaps)
	selectedSites []bool  // true for selected sites
	removegaps    bool    // If true, we will remove posision with >=1 gaps
	gamma         bool
	alpha         float64
}

func NewK2PModel(removegaps bool) *K2PModel {
	return &K2PModel{
		0,
		nil,
		removegaps,
		false,
		0.,
	}
}

/* computes K2P distance between 2 sequences */
func (m *K2PModel) Distance(seq1 []rune, seq2 []rune, weights []float64) (float64, error) {
	var dist float64

	trS, trV, _, _, total := countMutations(seq1, seq2, m.selectedSites, weights)
	trS, trV = trS/total, trV/total

	if m.gamma {
		dist = m.alpha * (.5*math.Pow(1.-2.*trS-trV, -1./m.alpha) + .25*math.Pow(1.-2.*trV, -1./m.alpha) - .75)
	} else {
		dist = -.5*math.Log(1.-2.*trS-trV) - .25*math.Log(1.-2.*trV)
	}

	return dist, nil
}

func (m *K2PModel) InitModel(al align.Alignment, weights []float64, gamma bool, alpha float64) (err error) {
	m.gamma = gamma
	m.alpha = alpha
	m.numSites, m.selectedSites = selectedSites(al, weights, m.removegaps)
	return
}
