package dna

import (
	"gonum.org/v1/gonum/mat"
)

type K2PModel struct {
	// Parameters (for eigen values/vectors computation)
	// Default 1.0
	// See https://en.wikipedia.org/wiki/Models_of_DNA_evolution#K80_model_(Kimura_1980)
	kappa float64
}

func NewK2PModel() *K2PModel {
	return &K2PModel{
		1.,
	}
}

// For Eigen values/vectors computation
//
func (m *K2PModel) InitModel(kappa float64) {
	m.kappa = kappa
}

func (m *K2PModel) Eigens() (val []float64, leftvectors, rightvectors *mat.Dense, err error) {
	val = []float64{
		0,
		-2 * (1 + m.kappa) / (m.kappa + 2),
		-2 * (1 + m.kappa) / (m.kappa + 2),
		-4 / (m.kappa + 2),
	}

	leftvectors = mat.NewDense(4, 4, []float64{
		1. / 4., 1. / 4., 1. / 4., 1. / 4.,
		0, 1. / 2., 0, -1. / 2.,
		1. / 2., 0, -1. / 2., 0,
		1. / 4., -1. / 4., 1. / 4., -1. / 4.,
	})

	rightvectors = mat.NewDense(4, 4, []float64{
		1., 0., 1., 1.,
		1., 1., 0., -1.,
		1., 0., -1., 1.,
		1., -1., 0., -1.,
	})

	return
}
