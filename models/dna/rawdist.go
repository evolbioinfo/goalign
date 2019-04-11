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
	countgapmut   bool    // If true, will count as 1 mutation '-' to 'A", default false
}

func NewRawDistModel(removegaps bool) *RawDistModel {
	return &RawDistModel{
		0,
		nil,
		removegaps,
		false,
	}
}

func (m *RawDistModel) SetCountGapMutations(countgapmut bool) {
	m.countgapmut = countgapmut
}

// computes the number of differences  between 2 sequences
// These differences include gaps vs. nt
func (m *RawDistModel) Distance(seq1 []rune, seq2 []rune, weights []float64) (diff float64, err error) {
	if m.countgapmut {
		diff, _ = countDiffsWithGaps(seq1, seq2, m.selectedSites, weights)
	} else {
		diff, _ = countDiffs(seq1, seq2, m.selectedSites, weights)
	}
	return
}

func (m *RawDistModel) InitModel(al align.Alignment, weights []float64) (err error) {
	m.numSites, m.selectedSites = selectedSites(al, weights, m.removegaps)
	return
}
