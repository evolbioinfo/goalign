package dna

import (
	"fmt"

	"github.com/evolbioinfo/goalign/align"
)

type PDistModel struct {
	numSites      float64 // Number of selected sites (no gaps)
	selectedSites []bool  // true for selected sites
	removegaps    bool    // If true, we will remove posision with >=1 gaps
	// If 0, will not count as 1 mutation '-' to 'A"
	// If 1, will count as 1 mutation '-' to 'A"
	// If 2, will count as 1 mutation '-' to 'A", but only the internal
	// Default 0
	countgapmut int
	// if true, we remove ambiguous positions for the normalisation by the length
	// default false
	// for example:
	// N vs. A : position not taken into account (can not decide wether there is a difference)
	// R vs. Y : position taken into account (we know there is a difference)
	removeAmbiguous bool
	sequenceCodes   [][]uint8 // Sequences converted into int codes
}

func NewPDistModel(removegaps bool) *PDistModel {
	return &PDistModel{
		0,
		nil,
		removegaps,
		0,
		false,
		nil,
	}
}

func (m *PDistModel) SetCountGapMutations(countgapmut int) (err error) {
	if countgapmut < 0 || countgapmut > 2 {
		err = fmt.Errorf("Gap count mode not available : %d", countgapmut)
	} else {
		m.countgapmut = countgapmut
	}
	return
}

// SetRemoveAmbiguous sets removeAmbiguous model variable
// if true, ambiguous positions are removed for the normalisation by the length
// for example:
// N vs. A : position not taken into account in length (can not decide wether there is a difference)
// R vs. Y : position taken into account in length (we know there is a difference)
func (m *PDistModel) SetRemoveAmbiguous(removeAmbiguous bool) {
	m.removeAmbiguous = removeAmbiguous
	return
}

/* computes p-distance between 2 sequences */
func (m *PDistModel) Distance(seq1 []uint8, seq2 []uint8, weights []float64) (diff float64, err error) {
	var total float64
	switch m.countgapmut {
	case GAP_COUNT_ALL:
		diff, total = countDiffsWithGaps(seq1, seq2, m.selectedSites, weights, m.removeAmbiguous)
	case GAP_COUNT_INTERNAL:
		diff, total = countDiffsWithInternalGaps(seq1, seq2, m.selectedSites, weights, m.removeAmbiguous)
	default:
		diff, total = countDiffs(seq1, seq2, m.selectedSites, weights, m.removeAmbiguous)
	}
	diff = diff / total
	return
}

func (m *PDistModel) InitModel(al align.Alignment, weights []float64, gamma bool, alpha float64) (err error) {
	m.numSites, m.selectedSites = selectedSites(al, weights, m.removegaps)
	m.sequenceCodes, err = alignmentToCodes(al)
	return
}

// Sequence returns the ith sequence of the alignment
// encoded in int
func (m *PDistModel) Sequence(i int) (seq []uint8, err error) {
	if i < 0 || i >= len(m.sequenceCodes) {
		err = fmt.Errorf("This sequence does not exist: %d", i)
		return
	}
	seq = m.sequenceCodes[i]
	return
}
