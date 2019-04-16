package dna

import (
	"gonum.org/v1/gonum/mat"
)

type JCModel struct {
}

func NewJCModel() *JCModel {
	return &JCModel{}
}

func (m *JCModel) InitModel() (err error) {
	return
}

// Left vectors and right vectors are given in column-major format
func (m *JCModel) Eigens() (val []float64, leftvectors, rightvectors *mat.Dense, err error) {
	val = []float64{
		0,
		-4. / 3.,
		-4. / 3.,
		-4. / 3.,
	}

	leftvectors = mat.NewDense(4, 4, []float64{
		1. / 4., 1. / 4., 1. / 4., 1. / 4.,
		-1. / 4., -1. / 4., 3. / 4., -1. / 4.,
		-1. / 4., 3. / 4., -1. / 4., -1. / 4.,
		3. / 4., -1. / 4., -1. / 4., -1. / 4.,
	})

	rightvectors = mat.NewDense(4, 4, []float64{
		1., 0., 0., 1.,
		1., 0., 1., 0.,
		1., 1., 0., 0.,
		1., -1., -1., -1.,
	})
	return
}
