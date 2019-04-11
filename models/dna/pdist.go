package dna

import (
	"github.com/evolbioinfo/goalign/align"
)

type PDistModel struct {
	numSites      float64 // Number of selected sites (no gaps)
	selectedSites []bool  // true for selected sites
	removegaps    bool    // If true, we will remove posision with >=1 gaps
	countgapmut   bool    // If true, will count as 1 mutation '-' to 'A", default false
}

func NewPDistModel(removegaps bool) *PDistModel {
	return &PDistModel{
		0,
		nil,
		removegaps,
		false,
	}
}

func (m *PDistModel) SetCountGapMutations(countgapmut bool) {
	m.countgapmut = countgapmut
}

/* computes p-distance between 2 sequences */
func (m *PDistModel) Distance(seq1 []rune, seq2 []rune, weights []float64) (diff float64, err error) {
	var total float64
	if m.countgapmut {
		diff, total = countDiffsWithGaps(seq1, seq2, m.selectedSites, weights)
	} else {
		diff, total = countDiffs(seq1, seq2, m.selectedSites, weights)
	}
	diff = diff / total
	return
}

func (m *PDistModel) InitModel(al align.Alignment, weights []float64) (err error) {
	m.numSites, m.selectedSites = selectedSites(al, weights, m.removegaps)
	return
}
