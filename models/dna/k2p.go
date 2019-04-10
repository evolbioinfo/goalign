package dna

import (
	"github.com/evolbioinfo/goalign/align"
	"math"
)

type K2PModel struct {
	numSites      float64 // Number of selected sites (no gaps)
	selectedSites []bool  // true for selected sites
	removegaps    bool    // If true, we will remove posision with >=1 gaps
	// Parameters (for eigen values/vectors computation)
	// Default 1.0
	// See https://en.wikipedia.org/wiki/Models_of_DNA_evolution#K80_model_(Kimura_1980)
	kappa float64
}

func NewK2PModel(removegaps bool) *K2PModel {
	return &K2PModel{
		0,
		nil,
		removegaps,
		1.,
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
		return 0, nil
	}
}

func (m *K2PModel) InitModel(al align.Alignment, weights []float64) (err error) {
	m.numSites, m.selectedSites = selectedSites(al, weights, m.removegaps)
	return
}

// For Eigen values/vectors computation
//
func (m *K2PModel) SetParameters(kappa float64) {
	m.kappa = kappa
}

func (m *K2PModel) Eigens() (val []float64, leftvector, rightvector [][]float64, err error) {
	val = []float64{
		0,
		-2 * (1 + m.kappa) / (m.kappa + 2),
		-2 * (1 + m.kappa) / (m.kappa + 2),
		-4 / (m.kappa + 2),
	}

	leftvector = [][]float64{
		[]float64{1. / 4., 1. / 4., 1. / 4., 1. / 4.},
		[]float64{0, 1. / 2., 0, -1. / 2.},
		[]float64{1. / 2., 0, -1. / 2., 0},
		[]float64{1. / 4., -1. / 4., 1. / 4., -1. / 4.},
	}

	rightvector = [][]float64{
		[]float64{1., 0., 1., 1.},
		[]float64{1., 1., 0., -1.},
		[]float64{1., 0., -1., 1.},
		[]float64{1., -1., 0., -1.},
	}
	return
}
