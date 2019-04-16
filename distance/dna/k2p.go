package dna

import (
	"math"

	"github.com/evolbioinfo/goalign/align"
)

type K2PModel struct {
	numSites      float64 // Number of selected sites (no gaps)
	selectedSites []bool  // true for selected sites
	removegaps    bool    // If true, we will remove posision with >=1 gaps
}

func NewK2PModel(removegaps bool) *K2PModel {
	return &K2PModel{
		0,
		nil,
		removegaps,
	}
}

/* computes K2P distance between 2 sequences */
func (m *K2PModel) Distance(seq1 []rune, seq2 []rune, weights []float64) (float64, error) {
	trS, trV, _, _, total := countMutations(seq1, seq2, m.selectedSites, weights)
	trS, trV = trS/total, trV/total

	dist := -0.5*math.Log(1.0-2.0*trS-trV) - 0.25*math.Log(1.0-2.0*trV)

	if dist > 0 {
		return dist, nil
	} else {
		return NT_DIST_MAX, nil
	}
}

func (m *K2PModel) InitModel(al align.Alignment, weights []float64) (err error) {
	m.numSites, m.selectedSites = selectedSites(al, weights, m.removegaps)
	return
}
