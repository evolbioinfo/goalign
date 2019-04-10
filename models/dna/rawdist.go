package dna

import (
	"github.com/evolbioinfo/goalign/align"
)

// Like pdist, but without
// Normalization by the number
// of sites
type RawDistModel struct {
	numSites      float64 // Number of selected sites (no gaps)
	selectedSites []bool  // true for selected sites
	removegaps    bool    // If true, we will remove posision with >=1 gaps
}

func NewRawDistModel(removegaps bool) *RawDistModel {
	return &RawDistModel{
		0,
		nil,
		removegaps,
	}
}

/* computes p-distance between 2 sequences */
func (m *RawDistModel) Distance(seq1 []rune, seq2 []rune, weights []float64) (float64, error) {
	diff, _ := countDiffs(seq1, seq2, m.selectedSites, weights)
	return diff, nil
}

func (m *RawDistModel) InitModel(al align.Alignment, weights []float64) (err error) {
	m.numSites, m.selectedSites = selectedSites(al, weights, m.removegaps)
	return
}
