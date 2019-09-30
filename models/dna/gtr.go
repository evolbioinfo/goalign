package dna

import (
	"fmt"

	"gonum.org/v1/gonum/mat"
)

type GTRModel struct {
	qmatrix    *mat.Dense
	leigenvect *mat.Dense
	val        []float64
	reigenvect *mat.Dense
}

func NewGTRModel() *GTRModel {
	return &GTRModel{
		nil,
		nil,
		nil,
		nil,
	}
}

//  /          \
// | *  d  f  b |
// | d  *  e  a |
// | f  e  *  c |
// | b  a  c  * |
//  \          /
func (m *GTRModel) InitModel(d, f, b, e, a, c, piA, piC, piG, piT float64) (err error) {
	m.qmatrix = mat.NewDense(4, 4, []float64{
		-(d*piC + f*piG + b*piT), d * piC, f * piG, b * piT,
		d * piA, -(d*piA + e*piG + a*piT), e * piG, a * piT,
		f * piA, e * piC, -(f*piA + e*piC + c*piT), c * piT,
		b * piA, a * piC, c * piG, -(b*piA + a*piC + c*piG),
	})
	// Normalization of Q
	norm := -piA*m.qmatrix.At(0, 0) -
		piC*m.qmatrix.At(1, 1) -
		piG*m.qmatrix.At(2, 2) -
		piT*m.qmatrix.At(3, 3)
	m.qmatrix.Apply(func(i, j int, v float64) float64 { return v / norm }, m.qmatrix)
	err = m.computeEigens()

	return
}

func (m *GTRModel) computeEigens() (err error) {
	var u mat.CDense

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
	eigen.VectorsTo(&u)
	reigenvect := mat.NewDense(4, 4, nil)
	leigenvect := mat.NewDense(4, 4, nil)
	reigenvect.Apply(func(i, j int, val float64) float64 { return real(u.At(i, j)) }, reigenvect)
	leigenvect.Inverse(reigenvect)

	m.leigenvect = leigenvect
	m.reigenvect = reigenvect
	m.val = val

	return
}

func (m *GTRModel) Eigens() (val []float64, leftvectors, rightvectors *mat.Dense, err error) {
	leftvectors = m.leigenvect
	rightvectors = m.reigenvect
	val = m.val

	return
}

func (m *GTRModel) Pij(i, j int, l float64) float64 {
	return -1
}

func (m *GTRModel) Analytical() bool {
	return false
}

func (m *GTRModel) NState() int {
	return 4
}
