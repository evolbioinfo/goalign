package dna

import (
	"fmt"

	"gonum.org/v1/gonum/mat"
)

type F81Model struct {
	// Parameters (for eigen values/vectors computation)
	// See https://en.wikipedia.org/wiki/Models_of_DNA_evolution#F81_model_(Felsenstein_1981)
	qmatrix    *mat.Dense
	leigenvect *mat.Dense
	val        []float64
	reigenvect *mat.Dense
}

func NewF81Model() *F81Model {
	return &F81Model{
		nil,
		nil,
		nil,
		nil,
	}
}

func (m *F81Model) InitModel(piA, piC, piG, piT float64) (err error) {
	m.qmatrix = mat.NewDense(4, 4, []float64{
		-(piC + piG + piT), piC, piG, piT,
		piA, -(piA + piG + piT), piG, piT,
		piA, piC, -(piA + piC + piT), piT,
		piA, piC, piG, -(piA + piC + piG),
	})
	// Normalization of Q
	norm := -piA*m.qmatrix.At(0, 0) -
		piC*m.qmatrix.At(1, 1) -
		piG*m.qmatrix.At(2, 2) -
		piT*m.qmatrix.At(3, 3)
	//norm := 1. / (2 * (piA*piC + piA*piG + piA*piT + piC*piG + piC*piT + piG*piT))
	m.qmatrix.Apply(func(i, j int, v float64) float64 { return v / norm }, m.qmatrix)

	//fmt.Printf("Q=%v\n", mat.Formatted(m.qmatrix, mat.Prefix("  "), mat.Squeeze()))

	err = m.computeEigens()

	return
}

func (m *F81Model) computeEigens() (err error) {
	// Compute eigen values, left and right eigenvectors of Q
	eigen := &mat.Eigen{}
	if ok := eigen.Factorize(m.qmatrix, mat.EigenRight); !ok {
		err = fmt.Errorf("Problem during matrix decomposition")
		return
	}

	val := make([]float64, 4)
	for i, b := range eigen.Values(nil) {
		val[i] = real(b)
	}
	u := eigen.VectorsTo(nil)
	reigenvect := mat.NewDense(4, 4, nil)
	leigenvect := mat.NewDense(4, 4, nil)
	reigenvect.Apply(func(i, j int, val float64) float64 { return real(u.At(i, j)) }, reigenvect)
	leigenvect.Inverse(reigenvect)

	m.leigenvect = leigenvect
	m.reigenvect = reigenvect
	m.val = val

	return
}

func (m *F81Model) Eigens() (val []float64, leftvectors, rightvectors *mat.Dense, err error) {
	leftvectors = m.leigenvect
	rightvectors = m.reigenvect
	val = m.val

	return
}

func (m *F81Model) pij(i, j int, l float64) float64 {
	return -1
}

func (m *F81Model) analytical() bool {
	return false
}
